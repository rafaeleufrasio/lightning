Function Header Example
=======================

.. highlight:: IDL

This is an example of a function header used in Lightning. In the example, all the details
needed to generate a function or procedure header for a new function are given. Here's what
the example would look like in the IDL source code, embedded in an example function it describes::

    function func, arg1, arg2, arg3, optional_input1=optional_input1, optional_input2=optional_input2, $
                   flag1=flag1, flag2=flag2, optional_output1=optional_output1, $
                   optional_output2=optional_output2
    ;+
    ; Name
    ; ----
    ;   FUNC
    ;
    ; Calling Sequence
    ; ----------------
    ;   ::
    ;
    ;       result = func(arg1, arg2, arg3 [, optional_input1 = , optional_input2 = , $
    ;                     optional_output1=optional_output1, optional_output2=optional_output2, $
    ;                     /flag1, /flag2])
    ;
    ; Inputs
    ; ------
    ;   arg1 : data type and size (e.g.: int array(2,2))
    ;       Note the space between the arg name, the colon, and the data type.
    ;       These lines should describe what the parameter is. It looks nice and
    ;       renders nicely if there's a hanging indent on this block. The lines
    ;       in this block must line up to render properly.
    ;   arg2 : data type and size (e.g.: float or double scalar)
    ;       This block should describe ``arg2``. References to variables or literal values
    ;       (i.e., ``3.14``) should be wrapped in double backquotes so they render in monospace
    ;       font.
    ;   arg3 : data type and size (e.g.: float or double array(N))
    ;       Description of ``arg3``. Note that inputs (and outputs) can have variable sizes.
    ;       This is indicated in the size description with a new variable. For example, 
    ;       ``arg3`` can be any 1-D array of length ``N``. If the number of array dimensions
    ;       is variable then an ellipsis (...) can be used (e.g., float or double array(...))
    ;
    ; Optional Inputs
    ; ---------------
    ;   optional_input1 : data type and size (e.g.: string array(M))
    ;       Description of ``optional_input1``. All optional inputs require a default value
    ;       at the end of the description outside of the final sentence period. The value itself 
    ;       should be wrapped in double backquotes as it is a literal value. This is the final
    ;       sentence at the end of which there is an example. (Default = ``[Hello, world]``)
    ;   optional_input2 : data type and size (e.g.: int scalar)
    ;       Description of ``optional_input2``. (Default = ``3``)
    ;   flag1 : data type and size (e.g.: flag)
    ;       Description of ``flag1``. Flags are simple parameters that can be set with ``/flag1``.
    ;       Setting the flag results in the parameter being equal to ``1``. Note that when not set,
    ;       the flag value is ``0``. Therefore, no default value is needed as by default they ``0``.
    ;   flag2 : data type and size (e.g.: flag)
    ;       Description of ``flag2``. Notice that the flag data type and size is simply flag.
    ;       Indicating that it is flag parameter.
    ;
    ; Output
    ; ------
    ;   result : data type and size (e.g.: int scalar)
    ;       This function outputs 0, because it isn't a real function. But in general here's where 
    ;       you'd describe what your function returns. If it is a procedure, which do not have a
    ;       dedicated returned variable, then the {variable name : data type and size} line above
    ;       can be removed, and the hanging indent on this description block can be removed too.
    ;
    ; Optional Outputs
    ; ----------------
    ;   optional_output1 : data type and size (e.g.: string scalar)
    ;       Description of ``optional_output1``.
    ;   optional_output2 : data type and size (e.g.: structure)
    ;       Description of ``optional_output2``. Notice that ``optional_output2`` is a
    ;       structure. In this case, it may be helpful to include a table of what the
    ;       structure contains. This can similarly be done for inputs and the main
    ;       output. Here is an example (and format) of what the structure description
    ;       table should look like. Note the the blank line between the table and this
    ;       description. Also, each column of the table is separated by five spaces.
    ;
    ;       ====     ============     ==============================================================================================
    ;       TAG      TYPE             DESCRIPTION
    ;       ====     ============     ==============================================================================================
    ;       TAG1     double(N, M)     First tag description. The tag column contains the tag name. 
    ;       TAG2     float            Second tag description. The type column has the format of type(size).
    ;       TAG3     string(5)        Third tag description. If the size in the type column is a scalar, then just include the type.
    ;       TAG4     string           Fourth tag description.
    ;       ====     ============     ==============================================================================================
    ;
    ; Notes
    ; -----
    ;   If a function had any additional notes that were not included in any other description,
    ;   here is where they would be. If any section is blank, for example there are no optional
    ;   outputs or notes, then the whole section can be removed from the header.
    ;
    ; Modification History
    ; --------------------
    ;   - 2022/03/17: Created (Erik B. Monson)
    ;   - 2022/06/23: Updated with final format (Erik B. Monson)
    ;   - YYYY/MM/DD: Something else (Somebody else)
    ;-

        return, 0

    end

And here is the rendered version of the example:

Name
----
  FUNC

Calling Sequence
----------------
  ::

      result = func(arg1, arg2, arg3 [, optional_input1 = , optional_input2 = , $
                    optional_output1=optional_output1, optional_output2=optional_output2, $
                    /flag1, /flag2])

Inputs
------
  arg1 : data type and size (e.g.: int array(2,2))
      Note the space between the arg name, the colon, and the data type.
      These lines should describe what the parameter is. It looks nice and
      renders nicely if there's a hanging indent on this block. The lines
      in this block must line up to render properly.
  arg2 : data type and size (e.g.: float or double scalar)
      This block should describe ``arg2``. References to variables or literal values
      (i.e., ``3.14``) should be wrapped in double backquotes so they render in monospace
      font.
  arg3 : data type and size (e.g.: float or double array(N))
      Description of ``arg3``. Note that inputs (and outputs) can have variable sizes.
      This is indicated in the size description with a new variable. For example, 
      ``arg3`` can be any 1-D array of length ``N``. If the number of array dimensions
      is variable then an ellipsis (...) can be used (e.g., float or double array(...))

Optional Inputs
---------------
  optional_input1 : data type and size (e.g.: string array(M))
      Description of ``optional_input1``. All optional inputs require a default value
      at the end of the description outside of the final sentence period. The value itself 
      should be wrapped in double backquotes as it is a literal value. This is the final
      sentence at the end of which there is an example. (Default = ``[Hello, world]``)
  optional_input2 : data type and size (e.g.: int scalar)
      Description of ``optional_input2``. (Default = ``3``)
  flag1 : data type and size (e.g.: flag)
      Description of ``flag1``. Flags are simple parameters that can be set with ``/flag1``.
      Setting the flag results in the parameter being equal to ``1``. Note that when not set,
      the flag value is ``0``. Therefore, no default value is needed as by default they ``0``.
  flag2 : data type and size (e.g.: flag)
      Description of ``flag2``. Notice that the flag data type and size is simply flag.
      Indicating that it is flag parameter.

Output
------
  result : data type and size (e.g.: int scalar)
      This function outputs 0, because it isn't a real function. But in general here's where 
      you'd describe what your function returns. If it is a procedure, which do not have a
      dedicated returned variable, then the {variable name : data type and size} line above
      can be removed, and the hanging indent on this description block can be removed too.

Optional Outputs
----------------
  optional_output1 : data type and size (e.g.: string scalar)
      Description of ``optional_output1``.
  optional_output2 : data type and size (e.g.: structure)
      Description of ``optional_output2``. Notice that ``optional_output2`` is a
      structure. In this case, it may be helpful to include a table of what the
      structure contains. This can similarly be done for inputs and the main
      output. Here is an example (and format) of what the structure description
      table should look like. Note the the blank line between the table and this
      description. Also, each column of the table is separated by five spaces.

      ====     ============     ==============================================================================================
      TAG      TYPE             DESCRIPTION
      ====     ============     ==============================================================================================
      TAG1     double(N, M)     First tag description. The tag column contains the tag name. 
      TAG2     float            Second tag description. The type column has the format of type(size).
      TAG3     string(5)        Third tag description. If the size in the type column is a scalar, then just include the type.
      TAG4     string           Fourth tag description.
      ====     ============     ==============================================================================================

Notes
-----
  If a function had any additional notes that were not included in any other description,
  here is where they would be. If any section is blank, for example there are no optional
  outputs or notes, then the whole section can be removed from the header.

Modification History
--------------------
  - 2022/03/17: Created (Erik B. Monson)
  - 2022/06/23: Updated with final format (Erik B. Monson)
  - YYYY/MM/DD: Something else (Somebody else)


Some additional notes:

This is kind of a hybrid between the IDL recommended header format and the python docs, especially numpy and astropy.

All the sections must be separated by a blank line to render correctly.

The example calling sequence used a multi-line literal block (i.e. things that you want to look like code, 
like long calling sequences) given by the specific syntax of the double colon. For smaller calling
sequences that are limited to a single line two backquotes can be placed around the calling sequence instead.
For example::

      ; Calling Sequence
      ; ----------------
      ;   ``my_procedure, arg1, arg2, flag1=flag1, optional_output1=optional_output1``
      ;

The above will render as:

Calling Sequence
----------------
  ``my_procedure, arg1, arg2, flag1=flag1, optional_output1=optional_output1``
