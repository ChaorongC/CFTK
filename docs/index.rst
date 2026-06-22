Cell-free DNA Toolkit Documentation
===================================
.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Start

   user_guide/model_power_calculator
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

CFTK is a versatile cfDNA analysis toolkit designed for cfDNA Bisulfite-sequencing data processing and analysing, aimed to support the biomarker discovery among lagre-scale Liquid Biopsy cohort. 

.. raw:: html

   <div style="text-align: center; margin: 20px 0;">
       <img src="_static/cftk_diagram.svg" alt="Power Analysis Overview" style="width: 700px; max-width: 100%;">
   </div>

We provide a model power calculator to evaluate whether a proposed biomarker discovery cohort is likely to produce a useful and detectable internally cross-validated classifier. 

.. grid:: 1
   :gutter: 2

   .. grid-item-card:: Model Power Calculator

      CFTK model power calculator, a tool to support the study design of cfDNA cohort.

      .. raw:: html

         <hr>
         <a href="https://cftk-model-power.streamlit.app" target="_blank" style="display: inline-block; background-color: #1b1233; color: white; padding: 10px 20px; text-decoration: none; border-radius: 4px; font-weight: bold; margin-bottom: 8px;">Try Now</a>
         <br>
         <a href="user_guide/model_power_calculator.html" style="color: #555; text-decoration: underline; font-size: 0.9em;">Tutorial &gt;</a>



We organized the CFTK workflow from rawdata processing to multimodal machine learning modeling and provide a report for results summarizing. The current workflow is driven by the ``cftk`` command implemented in ``cftk.py`` and controled by the user-prepared config file ``cftk_init.json``.

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
      :link: user_guide/index
      :link-type: doc

      Set up the initial config file ``cftk_init.json`` for your samples and explore standard CFTK workflow.

   .. grid-item-card:: Run Workflows
      :link: getting_started
      :link-type: doc

      Run CFTK workflow step by step or use the ``run-all`` command to finish all steps end to end.

   .. grid-item-card:: Command Reference
      :link: reference/cli
      :link-type: doc

      Explore all the available ``cftk`` commands and its function.


.. grid:: 1
   :gutter: 2

   .. grid-item-card:: Report Demo
      
      This report was generated using all the results from the CFTK workflow.
      
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
