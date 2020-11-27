Geometry
========

Basic types
-----------


.. class:: scine_utilities.ElementType

   Enum to represent element types including isotopes

   >>> h = ElementType.H # Represents isotopic mixture
   >>> h1 = ElementType.H1 # Represents only H with A = 1
   >>> d = ElementType.D # Represents only H with A = 2 (Deuterium)

.. autoclass:: scine_utilities.ElementInfo
   :members:

.. autoclass:: scine_utilities.Atom
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: scine_utilities.AtomCollection
   :members:
   :undoc-members:

   .. automethod:: __init__


Helper functions
----------------

.. automodule:: scine_utilities.geometry
   :members:
   :undoc-members:
