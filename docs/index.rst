Cell-free DNA Toolkit Documentation
===================================
.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Start

   installation
   getting_started

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: User Guide

   user_guide/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference

   reference/index

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Project

   development

.. rst-class:: cftk-hero

CFTK is a cfDNA analysis toolkit for bisulfite-sequencing processing and
downstream analysis, designed to support biomarker discovery in large
liquid-biopsy cohorts.

.. raw:: html

   <div style="text-align: center; margin: 20px 0;">
       <img src="_static/cftk_diagram.svg" alt="Power Analysis Overview" style="width: 700px; max-width: 100%;">
   </div>

We provide a model power calculator to evaluate whether a proposed biomarker discovery cohort is likely to produce a useful and detectable internally cross-validated classifier. 

.. grid:: 1
   :gutter: 2

   .. grid-item-card:: Model Power Calculator

      Run the model-development power calculator locally to support cfDNA cohort
      study design. No public hosted deployment is currently advertised.

      .. raw:: html

         <hr>
         <a href="user_guide/model_power_calculator.html" style="display: inline-block; background-color: #1b1233; color: white; padding: 10px 20px; text-decoration: none; border-radius: 4px; font-weight: bold; margin-bottom: 8px;">Run locally</a>

The beginner workflow starts with ``cftk init`` and ``cftk run``. CFTK validates
the schema-v2 project, runs core processing and QC with fail-fast stage
boundaries, and records exact commands, expected outputs, figures, tool
versions, and resume decisions. Advanced downstream analysis remains available
through explicit expert commands.

.. image:: _static/cftk_workflow.png
   :alt: CFTK Workflow Overview
   :align: center
   :width: 900px

Please follow the guides below to explore more details about the CFTK package.


.. grid:: 1 1 2 2
   :gutter: 2

   .. grid-item-card:: Installation
      :link: installation
      :link-type: doc

      Install CFTK, set up the enviroment and dependiencies for running. 

      
   .. grid-item-card:: Get Started
      :link: getting_started
      :link-type: doc

      Create a project, prepare the default reference profile, and run the
      validated beginner workflow.

   .. grid-item-card:: Workflow Guides
      :link: user_guide/index
      :link-type: doc

      Inspect expected outputs and use expert workflows for processing, QC, and
      downstream analysis.

   .. grid-item-card:: Command Reference
      :link: reference/cli
      :link-type: doc

      Explore all the available ``cftk`` commands and its function.


.. grid:: 1
   :gutter: 2

   .. grid-item-card:: Report Demo
      
      Static legacy report preview for layout and navigation only. It is not a
      current default report or a patient-data result.
      
      .. raw:: html

         <hr>
         <a href="_static/sample_report.html" target="_blank" style="
            display: inline-block;
            background-color: #1b1233;
            color: white;
            padding: 10px 20px;
            text-decoration: none;
            border-radius: 4px;
            font-weight: bold;
            margin-bottom: 8px;
         ">Open report</a>
         <br>
         <a href="_static/sample_report.html" download style="color: #555; text-decoration: underline; font-size: 0.9em;">Download full report &gt;</a>
