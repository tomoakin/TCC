TCC
===


Title
-----

TCC: Differential expression analysis for tag count data with robust normalization strategies

Description
----------------------------------------------------------------------------------------

This package provides a series of functions for performing
differential expression analysis from RNA-seq count data using robust
normalization strategy (called DEGES). The basic idea of DEGES is that
potential differentially expressed genes or transcripts (DEGs) among
compared samples should be removed before data normalization to obtain
a well-ranked gene list where true DEGs are top-ranked and non-DEGs are
bottom ranked. This can be done by performing a multi-step normalization
strategy (called DEGES for DEG elimination strategy). A major
characteristic of TCC is to provide the robust normalization methods for
several kinds of count data (two-group with or without replicates,
multi-group/multi-factor, and so on) by virtue of the use of combinations
of functions in other sophisticated packages (especially edgeR, DESeq,
and baySeq).

Documentation
-------------

Documents of this software are available at
http://www.iu.a.u-tokyo.ac.jp/~kadota/TCC/  
in addition to the documents included in the software package.

Development versions
--------------------

The development version is now on the github repository
at https://github.com/kohijiri/huaying.asagao

Releases
--------

Release of this software is made through CRAN
at http://cran.r-project.org/web/packages/TCC/


No warranty
-----------

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

Lisense
-------

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License 
version 2 (GPL-2).

A copy of GPL-2 is accompanied in this directory in a file named COPYING.

Authors
-------

Jianqiang Sun, Tomoaki Nishiyama, Kentaro Shimizu, and Koji Kadota
