Multi-modal Modeling
====================================

MESA Modeling
-------------

CFTK includes MESA-style multimodal modeling commands.

Run modality performance screening:

.. code-block:: bash

   cftk --config cftk_init.json mesa --performance

Run model construction:

.. code-block:: bash

   cftk --config cftk_init.json mesa --mesa-model

Run leave-one-out cross-validation and plots:

.. code-block:: bash

   cftk --config cftk_init.json mesa --loocv

Commonly used together:

.. code-block:: bash

   cftk --config cftk_init.json mesa --performance --mesa-model --loocv
