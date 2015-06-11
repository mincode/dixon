Installation instructions for
# dr.mpl
# Package for Dixon Resultant Computation 
# (C) 2015 Manfred Minimair
# Licensed under the GPL 3.0 (http://www.gnu.org), see gpl-3.0.txt distributed with this package

Import dr.mpl into your Maple session using the command read "dr.mpl";
The Dixon resultant function will be made accessible as
DR:-DixonResultant(...)

Usage:
  DR:-DixonResultant(plist, vlist)
  # plist # list of inhomogeneous polynomials in the variables in vlist
  # vlist # list of variables, nops(plist) = nops(vlist) + 1
  # Returns the Dixon resultant of the polynomials in plist.

Open dr.mws for further help on the package.

How to cite the package:
BibTeX entry:
@Misc{minimair15:_dr,
  author =	 {M. Minimair},
  title =	 {DR: Dixon Resultant Package for Maple},
  howpublished = {http://minimair.org/dr},
  year =	 2015
}

