DecompFlux
==========

The m-files are for decomposing flux distributions into sets of elementary flux modes (EFMs) without finding the EFM matrix in prior. 'decompflux.m' is the main function. The principle of the algorithm has been described in Chan and Ji (2011).

Contents:
A.	Functions
B.	Optimization solver
C.	Technical problems
D.	Example
E.	Reference

A. Functions
1.	'decompflux.m' is the main function file. All other functions are invoked by it. Given the stoichiometrix matrix, a flux distribution can be decomposed into a set of EFM by this function. Parameters have been default if not supplied by user.
If you simply want an arbitrary decomposition, you can set CT=2 and do not give any objective function. If you want EFMs with possibly larger contributions to a particular reaction, you can try both CT=1 and 2. In general, CT=1 should give better results but you may have to adjust the parameter (see section ‘Technical problems’ point 3).
2.	‘isEFM.m’ is a function to check whether a set of flux modes is EFMs or not given the stoichiometric matrix by solving LP. It invokes ‘isDC.m’. For flux distributions which can cause numerical difficulties, say, when involving non-integer, and especially when the magnitude differs largely between fluxes of different reactions (e.g. the largest flux is 5000 and the smallest flux is 0.1), the practical solution of the MILP by a solver may be prone to error because of tolerance problems. Checking by LP is thus a good way to make sure all the modes are elementary.
3.	‘SreToSir.m’ is a simple function to transform a stoichiometrix matrix with reversible reactions into one without, by splitting each reversible reaction into two irreversible reactions.
**'decompflux.m' requires an MILP solver. The current version supports Cplex11 or above or solvers embedded in the COBRA toolbox if the COBRA toolbox has been installed. ‘isEFM.m’ requires an LP solver. The current version supports Cplex11 or above or solvers embedded in the COBRA toolbox or the ‘linprog.m’ native in Matlab. Please see the section ‘Optimization solver’ below for further information.

B. Optimization solver
The algorithms in ‘decompflux.m’ and ‘isEFM.m’ are written so that they run within Matlab with the aid of Cplex 11 or above or any solvers called by the COBRA toolbox. You can change the code to suit your own optimization software. In this case, you should change the m file 'DFD1sub.m' etc. as 'DFD1.m' includes only the processing of flux modes, but not optimization.
•	Cplex:
If you want to use Cplex and do not have one, it is very simple. IBM offers Cplex for free for people in academic field. You can go to join the IBM Academic Initiative (or just google the term):
https://www.ibm.com/developerworks/university/membership/join.html
After joining, you can go to:
https://www14.software.ibm.com/webapp/iwm/web/reg/signup.do?source=scholars
to log in and search for Cplex to download.
Once finished, you just add the path containing the Cplex function for Matlab, then it can work.
(e.g. C:\Program Files\IBM\ILOG\CPLEX_Studio124\cplex\matlab\x64_win64)
•	COBRA toolbox:
COBRA toolbox is a powerful toolbox integrated with various computational methods based on constriant-based modeling of metabolic network. It is available at:
http://opencobra.sourceforge.net
You need to install it and attach it to Matlab. Also, you need to attach solver to the toolbox for it to work. Gurobi (http://www.gurobi.com/) is a powerful optimization software free for academic purpose and supported by COBRA toolbox.
The attached ‘solveCobraMILPedit.m’ is one of the m-file in the COBRA toolbox edited by the author because the original one contains error when using Gurobi 5 as the solver. The script is part of the openCOBRA Project and is distributed under the GNU GPLv3 or later.  However, this software is designed for scientific research and as such may contain algorithms that are associated with patents in the U.S. and abroad.  If the user so chooses to use the software provided by the openCOBRA project for commercial endeavors then it is solely the user’s responsibility to license 
any patents that may exist and respond in full to any legal actions taken by the patent holder.
For details, please visit http://opencobra.sourceforge.net/

C. Technical problems
If you have problems in obtaining a correct decomposition, here are some possible reasons:
1.	Any reversible reactions in the stoichiometrix matrix/ flux distribution? Transform them into irreversible one.
2.	If you have input an objective function, have you made sure the sign is correct? Minimizing uses ‘+’ and maximizing uses ‘-‘. Also, since each reversible reaction appears as two individual reactions, the sum of all fluxes should be minimized to get meaningful results. Thus, an example objective vector is ‘-10000’ for biomass reaction (max biomass) and ‘1’ for all others (min sum of fluxes).
3.	Get non-EFM flux modes in the solution? 
Realistic or genome-scale flux distributions can cause numerical difficulties, for example, when involving non-integer, and especially when the magnitude differs largely between fluxes of different reactions (e.g. the largest flux is 5000 and the smallest flux is 0.1). The practical solution of the MILP by a solver may be prone to error because of tolerance problems.
•	Try lowering the tolerance for integer feasibility, primal feasibility, dual feasibility (or optimality) and turn off the automatic scaling of the model.
•	Try to scale the flux distribution so that the smallest non-zero flux is in the order of 0.1.
(At the same time, the largest flux cannot be too large which also causes numerical problem. Flux distribution with fluxes within the range of 0.1 to 106 has been successfully decomposed.
•	Try adjusting the large number ‘big’ in the input and use ‘isEFM.m’ to check. Integer programming solver can sometimes be quite sensitive to that number. Using CT=2 in ‘decompflux.m’ in general has no problem. But when using CT=1 (to approximate the largest contributing EFM), especially if your flux distribution have fluxes very diverse in their magnitudes, problems may appear. You should try to vary the order of magnitude of the large number ‘big’.
•	Decompose the non-elementary modes again using ‘decompflux.m’.


D. Example
The script file 'test.m' can run an example contained in 'test.mat'. It is the flux distribution describing the growth of E. coli on the complex LB medium analyzed in Chan and Ji (2011). It is a hard example because the number of non-zero entries in the flux distribution is large and the magnitude of fluxes varies a lot with the largest being 70000 and the smallest being 0.1 after scaling. It should results at a decomposition of 25 EFMs as in Chan and Ji (2011).
The genome-scale metabolic network used here is the E. coli iAF1260 network originating from Feist (2007). ‘FD.S’, ‘FD.rxn’ and ‘FD.met’ are the stoichiometric matrix, reaction and metabolite names reproduced from the original network which is freely available in the BiGG database: http://bigg.ucsd.edu/
It is distributed according to the following license:
Copyright 2007 The Regents of the University of California
All Rights Reserved
Permission to use, copy, modify and distribute any part of this BiGG Database for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
Those desiring to incorporate this BiGG database into commercial products or use for commercial purposes should contact the Technology Transfer & Intellectual Property Services, University of California, San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910, Ph: (858) 534-5815, FAX: (858) 534-7345, e-mail:invent@ucsd.edu.
IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS BIGG DATABASE, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
THE BIGG DATABASE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. THE UNIVERSITYOF CALIFORNIA MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE BIGG DATABASE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.

E. Reference
Chan,S.H.J. and Ji,P. (2011) Decomposing flux distributions into elementary flux modes in genome-scale metabolic networks. Bioinformatics, 27(16), 2256-2262.
Feist,A.M. et al. (2007) A genome-scale metabolic reconstruction for Escherichia coli K-12 MG1655 that accounts for 1260 ORFs and thermodynamic information. Mol. Syst. Biol., 3, 121.

Please cite the paper if you has used the algorithm in your publication.
I earnestly hope it can help you. Please feel free to contact me if you have any question.
Thank you!

Siu Hung Joshua Chan
email: joshua.chan@connect.polyu.hk
