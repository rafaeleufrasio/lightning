==========================
Frequently Asked Questions
==========================
 
How to Easily Install Dependencies
----------------------------------

Lightning requires the `IDL Astronomy User's Library <https://idlastro.gsfc.nasa.gov>`_, `IDL Coyote <http://www.idlcoyote.com>`_, and
`Craig Markwardt's MPFIT library <http://purl.com/net/mpfit>`_. All three will need to be installed and `added to your IDL
path <https://www.l3harrisgeospatial.com/Support/Self-Help-Tools/Help-Articles/Help-Articles-Detail/ArtMID/10220/ArticleID/16156/Quick-tips-for-customizing-your-IDL-program-search-path>`_
before using Lightning. To easily download the dependencies, we recommend making a new directory, changing into the new
directory, and using the following shell commands::

    git clone https://github.com/wlandsman/IDLAstro.git
    git clone https://github.com/idl-coyote/coyote.git
    wget https://pages.physics.wisc.edu/~craigm/idl/down/mpfit.tar.gz
    mkdir mpfit
    tar -xf mpfit.tar.gz -C mpfit
    rm mpfit.tar.gz
