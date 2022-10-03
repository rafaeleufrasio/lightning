LIGHTNING_CONSTANTS
===================

Name
----
LIGHTNING_CONSTANTS

Purpose
-------
Generates a system variable that contains all of the constants needed
for unit conversions in Lightning. All constants are in CGS units.

Calling Sequence
----------------
::

    lightning_constants

Output
------
The created system variable, ``!lightning_cgs``.

Modification History
--------------------
- 2016/05/01: Created (Rafael T. Eufrasio)
- 2022/03/22: Changed file name and made into a procedure (Keith Doore)
- 2022/03/22: Added documentation (Keith Doore)
- 2022/04/18: Removed structure name to prevent multiple call issues (Keith Doore)
- 2022/05/09: Changed system variable name from ``!cv`` to ``!lightning_cgs`` (Keith Doore)

