Visualization And Reports
====================================

Visualization Generation
--------------------------
In standard CFTK workflow, visualization for the results were output automaticlly. You can also use ``vis`` to regenerate figures from existing result files:

.. code-block:: bash

   cftk --config cftk_init.json vis --mode all
   cftk --config cftk_init.json vis --mode power qc diff

Report Generation
-----------------

Generate a HTML report:

.. code-block:: bash

   cftk --config cftk_init.json report

Reports are written under:

.. code-block:: text

   <output_dir>/results/report/
