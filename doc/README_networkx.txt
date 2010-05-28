NOTE
====

In order for the CombinedGraph plotting to work one has to patch any version
of networkx < 0.37.

Apply the two patches:

  cd networkx
  patch -p0 < ...../edgecolor_OB20080124.patch
  patch -p0 < ...../nodelinewidths_OB20080124.patch

Newer versions of networkx incorporate these patches.

