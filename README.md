**Note: some functions may require additional licenses from QuantumBio, CCG, or other 3rd parties.**

# MOE/DivCon

The MOE/DivCon Scientific Vector Language (SVL) package is a set of support applications which 
allow QuantumBio software to be installed in and used with the MOE and MOE/batch platform from 
[Chemical Computing Group, Inc.](https://www.chemcomp.com/).

The software consists of a self-contained set of SVL applications and scripts and an associated menu 
file which can be installed in a user's home directory using the following steps.

A more detailed set of installation instructions (along with tutorials) can be found on the 
QuantumBio [WebSite](http://www.quantumbioinc.com/resources/manual/installation/).

## The ${HOME}/moefiles directory 

As detailed in the MOE manual, you may install software and plugins such as those provided by 
QuantumBio using the moefiles directory placed in the user's home directory. This is the easiest 
place to install these plugins as the user generally does not need admin privileges to do so.

```
    % mkdir -p ${HOME}/moefiles/svl
    % mkdir -p ${HOME}/moefiles/menu
```

## Checkout and Install the included SVL directory

You are welcome to clone the entirety of the MOE/DivCon source tree at a location of your own choosing. 
We will then create a symbolic link from within the ${HOME}/moefiles/svl to that location.

```
    % cd /path/you/wish/to/install/MOEDivCon
    % git clone --branch GridMarkets https://github.com/quantumbio/MOEDivCon.git MOEDivCon-GridMarkets
    % cd ${HOME}/moefiles/svl
    % ln -s /path/you/wish/to/install/MOEDivCon/MOEDivCon-GridMarkets/svl qb.svl
```

## Install the included QuantumBio Menu

Finally, to easily access the software from within MOE, you should install the QuantumBio menu. 
This menu can be found in the Extra menu on the main MOE window. 

**After you install MOE/DivCon using these steps, restart MOE and MOE/batch.**

```
    % cd ${HOME}/moefiles/menu
    % ln -s /path/you/wish/to/install/MOEDivCon/MOEDivCon-GridMarkets/menu/menu-extra-MTScore menu-extra-MTScore
```

![Extra Menu](./doc/images/qb_extra_menu.png)

## GridMarkets Envoy Installation

When you have completed the installation of the QuantumBio-GridMarkets Integration Manager, if you have not already done so, sign up for a GridMarkets account and install the Envoy package on your local workstation/laptop. This software will take care of the secure communications between your environment and the GridMarkets platform. The instructions for this process are found on the following link:

 * https://www.pharma.gridmarkets.com/quantumbio-setup

## Examples

When you have successfully installed our QuantumBio-GridMarkets Interface Manager and the associated Envoy tool, you can execute GridMarkets calculations remotely (on the cloud) just as you would run them locally. A MOE/DivCon based tutorial is available on the following link (and webinars will be added over comming weeks and months).

 * https://www.pharma.gridmarkets.com/quantumbio
 
The Manager also includes a command line based tool which is a modified version of the qmechanic (DivCon) wrapper script supplied in our regular distribution. This script supports the addition of the --cloud command line option which communicates with Envoy (which ultimately communicates with the GridMarkets cloud).
 
  * General QM/MM: http://www.quantumbioinc.com/resources/manual/divconcli/
  * MovableType: http://www.quantumbioinc.com/resources/manual/movabletype/

For example, this command line will run the MTScoreES (EndState) simulation on a single protein:ligand complex as referenced in the following online tutorial:

 * https://www.quantumbioinc.com/resources/manual/movabletype/#MTScoreES-cli

```
$ wget http://downloads.quantumbioinc.com/media/tutorials/MT/4w7t_protein.pdb
$ wget http://downloads.quantumbioinc.com/media/tutorials/MT/4w7t_ligand.mol2
$ /path/to/MOEDivCon-GridMarkets/bin/qmechanic 4w7t_protein.pdb --ligand 4w7t_ligand.mol2 --mtscore endstate --cloud GridMarkets
```

Any of the DivCon/qmechanic tutorials found on the QuantumBio website may be run with the addition of --cloud GridMarkets on the command line. Please see the simple pricing for this cloud-based service on this link:

 * https://www.pharma.gridmarkets.com/quantumbio
 
