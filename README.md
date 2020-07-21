# Benredjem-Gallion

This code was originally used in the publication:
Exploring use of unsupervised clustering to associate signaling profiles of GPCR ligands to clinical response Nat Commun . 2019 Sep 9;10(1):4075. doi: 10.1038/s41467-019-11875-6.

This publication used a series of .py scripts run in series. This code base has since been updated to jupyter notebooks, but achieves the same purpose.

Purpose: Utilize multiple iterations of NNMF followed by kmeans in order to define a similarity score for a series of compounds based on the dose response curves for a series of biological assays. NNMF reduces the dimensionality of the data into representative basis vectors, kmeans then uses these basis vectors in order to cluster. However, since both NNMF and Kmeans start off with a randomization, multiple iterations can converge to different minima. By running this process multiple times with multiple K values we can quantify the frequency with which two compounds cluster together, and therefore assign a similarity score in their original dose response. For a more detailed explaination please refer to the method of the paper above.

Input: A structured csv or xlsx file with unique compounds on each row and each column representing a dose response quantification (e.g. Emax, LogR, etc) across a series of assays. In the case of the original paper (example data included) 25 compounds targetting GPCRs were tested against a battery of 10 assays generating a dose response curve for each. This dose response curve was then characterized (using external software) into representative quantifications (e.g. Emax, pEC50, LogR, pLogKA, Log(T/Ka)).

Output: A compound x compound frequency matrix indicating the frequency that each compound coclustered across the entire method. These values can be interpretted as a similiarity score.


Pipepline Walkthough.

Step1
