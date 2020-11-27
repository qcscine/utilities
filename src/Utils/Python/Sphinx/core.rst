Core
====


The ``core`` submodule defines interfaces that permeate the SCINE Project.
Individual components may offer models of these interfaces that are then
available through the :class:`core.ModuleManager` class.

ModuleManager
-------------

.. autoclass:: scine_utilities.core.ModuleManager
   :members:
   :undoc-members:

Interfaces
----------

Calculator
~~~~~~~~~~

.. autoclass:: scine_utilities.core.Calculator
   :members:
   :undoc-members:

.. autofunction:: scine_utilities.core.get_available_settings

.. autofunction:: scine_utilities.core.get_possible_properties

.. autofunction:: scine_utilities.core.load_system

.. autofunction:: scine_utilities.core.load_system_into_calculator

Logging
-------

.. autoclass:: scine_utilities.core.Log
   :members:
   :undoc-members:

.. autoclass:: scine_utilities.core.Sink
