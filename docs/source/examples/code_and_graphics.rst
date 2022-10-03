Code and Graphics
=================

.. highlight:: idl

Here's what a ``doctest`` block looks like, for reference; these are intended for showing interactive use in
python documentation.

>>> 2 + 2
5

For IDL we should probably use the ``code-block`` directive which allows us to more explicitly specify
the syntax highlighting

.. code-block:: idl

    ; Comments get italicized and a different color.
    print, 2 + 2 ; Comments after code lines don't look right though

You can also use the default ``::`` syntax for creating code blocks, which requires you to specify the ``highlight``
direcive elsewhere in the document::

    ; Comments get italicized and a different color.
    print, 2 + 2 ; Comments after code lines don't look right though

If we want to fake a command prompt to make it clear that we're showing something interactive at the IDL prompt, we can
just type in a fake prompt, but it looks weird because the ``>`` is a different color, and comments won't work as I
showed above::

    IDL> print, 'The > looks funny.'
    IDL> ; and comment lines aren't highlighted properly (note that keywords "and" strings are highlighted)

To put in images, you can use the ``image`` directive, and specify the absolute path of the image relative to the
top of the documentation source directory:

.. image:: /_static/logo_black.png
    :align: center
    :width: 100px

There are the typical options you would expect, like alignment and width.

There's also a ``figure`` directive, which lets you have a caption and whatnot:

.. figure:: /_static/logo_black.png
    :align: center
    :width: 100px

    Here we show the logo for Lightning, which is the only image in the documentation
    directory. I don't *think* that you can label the figures as Figure 1, etc., but
    I might just have missed it. You can, however, use the reST ``target`` system to
    put a target in the figure caption and link to it from elsewhere on the page.
