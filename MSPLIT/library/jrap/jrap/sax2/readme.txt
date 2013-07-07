Overview
--------
JRAP is a quick port of Patrick Pedrioli's RAP program for
reading MSXML data files in a random access fashion.  It
has not been fully tested and serves more as a guide
to managing the Xerces XML parser on an indexed XML
file.

An example program which uses these classes has been 
provided "TestParser".  TestParser reads a MSXML file
and prints the contents of one or more scans which
you specify on the command line ie.:

  java TestParser msXMLfilename Scan# <Scan#> ...
  
The SpectrumViewer in the examples package is another small programs
that shows the use of the parser. It displays spectra from a 
mzXML file graphically. It needs some additional jar 
files, it explains how to get them by running java -jar jrap.jar


Files:
jakarta-regexp-1.4.jar - regular expression tools required for jrap_3.0.jar;
needed in your classpath (added 20070227)
jrap_3.0.jar - update of JRAP for mzXML 3.0 support (added 20070227)
jrap-dist.jar
