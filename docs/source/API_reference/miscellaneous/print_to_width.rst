PRINT_TO_WIDTH
==============

Name
----
PRINT_TO_WIDTH

Purpose
-------
Prints text to the terminal or IDL workbench. Strings that
are longer than the width of the window will be wrapped to
the next line without splitting words. This is different
from IDL's ``print`` procedure, which will split words.

Calling Sequence
----------------
::

    print_to_width, str

Input
-----
``str`` : string scalar
    The string that is to be printed.

Output
------
The input string is printed to the screen, without splitting words.

Modification History
--------------------
- 2022/05/12: Created (Keith Doore)

