#!/usr/bin/env bash

set -euo pipefail

if [[ "$#" -lt 4 || "$#" -gt 5 ]]; then
    echo "Usage: $0 MANIFEST REFERENCE_FASTA TARGET_INTERVAL_LIST OUTPUT_DIR [THREADS]" >&2
    exit 2
fi

manifest=$1
reference_fasta=$2
target_interval_list=$3
output_dir=$4
threads=${5:-10}

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "${script_dir}/../.." && pwd)
compare_script="${script_dir}/compare_duplicate_marking.py"

test -s "${manifest}"
test -s "${reference_fasta}"
test -s "${reference_fasta}.fai"
test -s "${target_interval_list}"
test -f "${compare_script}"
if [[ -e "${output_dir}" ]]; then
    echo "Refusing to overwrite existing validation output: ${output_dir}" >&2
    exit 2
fi

mkdir -p "${output_dir}" "${output_dir}/markdup" "${output_dir}/metrics" \
    "${output_dir}/methylation" "${output_dir}/resources" "${output_dir}/tmp"
cp "${manifest}" "${output_dir}/input_manifest.tsv"

exec 19>"${output_dir}/commands.trace"
export BASH_XTRACEFD=19
export PS4='+ ${BASH_SOURCE}:${LINENO}: '
set -x

{
    date --iso-8601=seconds
    python --version
    samtools --version | head -1
    sambamba --version | sed -n '2p'
    picard MarkDuplicates --version || test "$?" -eq 1
    MethylDackel --version
} >"${output_dir}/tool_versions.txt" 2>&1

comparison_manifest="${output_dir}/comparison_manifest.tsv"
printf 'sample\tgroup\tinput_bam\tsambamba_bam\tpicard_bam\tsambamba_metrics\n' \
    >"${comparison_manifest}"

tail -n +2 "${manifest}" | while IFS=$'\t' read -r sample group input_bam ot_bounds ob_bounds; do
    test -n "${sample}"
    test -s "${input_bam}"
    sample_markdup="${output_dir}/markdup/${sample}"
    sample_metrics="${output_dir}/metrics/${sample}"
    sample_methylation="${output_dir}/methylation/${sample}"
    sample_tmp="${output_dir}/tmp/${sample}"
    mkdir -p "${sample_markdup}" "${sample_metrics}" "${sample_methylation}" "${sample_tmp}"

    sambamba_bam="${sample_markdup}/sambamba.bam"
    picard_bam="${sample_markdup}/picard.bam"
    sambamba_markdup_metrics="${sample_markdup}/sambamba.metrics.txt"
    picard_markdup_metrics="${sample_markdup}/picard.metrics.txt"

    /usr/bin/time -v -o "${output_dir}/resources/${sample}.sambamba.time.txt" \
        sambamba markdup -t "${threads}" "${input_bam}" "${sambamba_bam}" \
        2>"${sambamba_markdup_metrics}"
    samtools index -@ "${threads}" "${sambamba_bam}"

    TMPDIR="${sample_tmp}" /usr/bin/time -v \
        -o "${output_dir}/resources/${sample}.picard.time.txt" \
        picard -Xmx8g MarkDuplicates \
        I="${input_bam}" O="${picard_bam}" R="${reference_fasta}" \
        M="${picard_markdup_metrics}" SORTING_COLLECTION_SIZE_RATIO=0.15 \
        ASSUME_SORT_ORDER=coordinate OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        MAX_RECORDS_IN_RAM=1000 TMP_DIR="${sample_tmp}"
    samtools index -@ "${threads}" "${picard_bam}"

    for method in sambamba picard; do
        method_bam="${sample_markdup}/${method}.bam"
        /usr/bin/time -v -o "${output_dir}/resources/${sample}.${method}.hs_metrics.time.txt" \
            picard -Xmx8g CollectHsMetrics \
            I="${method_bam}" O="${sample_metrics}/${method}.hs_metrics.txt" \
            R="${reference_fasta}" BAIT_INTERVALS="${target_interval_list}" \
            TARGET_INTERVALS="${target_interval_list}" \
            PER_TARGET_COVERAGE="${sample_metrics}/${method}.per_target_coverage.txt" \
            MINIMUM_MAPPING_QUALITY=20 COVERAGE_CAP=1000 NEAR_DISTANCE=500

        /usr/bin/time -v -o "${output_dir}/resources/${sample}.${method}.methylation.time.txt" \
            MethylDackel extract --mergeContext --minDepth 10 \
            --maxVariantFrac 0.25 -@ "${threads}" --OT "${ot_bounds}" \
            --OB "${ob_bounds}" -o "${sample_methylation}/${method}" \
            "${reference_fasta}" "${method_bam}"
    done

    samtools quickcheck -v "${sambamba_bam}" "${picard_bam}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${sample}" "${group}" "${input_bam}" "${sambamba_bam}" \
        "${picard_bam}" "${sambamba_markdup_metrics}" >>"${comparison_manifest}"
done

PYTHONPATH="${repo_root}/src:${repo_root}" python "${compare_script}" \
    "${comparison_manifest}" "${output_dir}/duplicate_marking_comparison.json" \
    --metrics-dir "${output_dir}/metrics" \
    --methylation-dir "${output_dir}/methylation"

set +x
exec 19>&-
find "${output_dir}" -type f ! -name checksums.sha256 -print0 \
    | sort -z | xargs -0 sha256sum >"${output_dir}/checksums.sha256"
echo "Duplicate-marking comparison completed: ${output_dir}"
