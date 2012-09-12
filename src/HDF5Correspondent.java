/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.ListIterator;
import java.util.Vector;
import javax.swing.tree.DefaultMutableTreeNode;
import ncsa.hdf.object.Datatype;
import ncsa.hdf.object.FileFormat;
import ncsa.hdf.object.HObject;
import ncsa.hdf.object.h5.H5CompoundDS;
import ncsa.hdf.object.h5.H5Group;
import ncsa.hdf.object.h5.H5ScalarDS;
import ncsa.hdf.object.h5.H5File;
import svljava.SVLJava;
import svljava.SVLJavaException;
import svljava.SVLVar;
import svljava.SVLWriter;

        
/**
 *
 * @author roger
 */
public class HDF5Correspondent implements SVLJavaDriver {

    private int taskID = 0;		// SVL task identifier
    
    @Override
    public SVLVar run(SVLVar arg, int tid) throws SVLJavaException
    {
	taskID = tid;

	try {
	    String cmd = arg.getTokn(1);	// subcommand
            SVLJava.print("cmd="+cmd);
	    SVLVar res = arg.peek(2);		// its argument
            SVLJava.print("length="+res.length());
	    if(cmd.equals("retrieveQMScore" ))
            {
                res = retrieveQMScore(res);
            }
	    else if (cmd.equals("retrieveResiduePWD"))
            {
                res = retrieveResiduePWD(res);
            }
	    else if (cmd.equals("retrieveAtomByAtomPWD"))
            {
                res = retrieveAtomByAtomPWD(res);
            }
	    else {
		throw new SVLJavaException("Unknown command: '" + cmd + "'.");
	    }
	    return res;
	}
	catch (Exception ex) {
	    StringWriter sw = new StringWriter();
	    ex.printStackTrace(new PrintWriter(sw));
	    throw new SVLJavaException("Exception:\n"+sw.toString());
	}
    }

    private SVLVar retrieveQMScore(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
        String xPath="/DivCon/"+target+"/QM Score";
        HObject ho=findHDF5Object(h5File, xPath);
        H5CompoundDS hObject=(H5CompoundDS)ho;
        if(hObject==null) throw new SVLJavaException("QM Score does not exist for: '" + target + "'.");
        hObject.init();
        hObject.clear();
        hObject.setMemberSelection(true);
        Vector o=(Vector)hObject.getData();
        String[] dimNames=hObject.getDimNames();
    SVLVar[] data=new SVLVar[hObject.getMemberCount()];
        for(int index=0;index<hObject.getMemberCount();index++)
        {
            hObject.selectMember(index);
            if(index==0)
            {
                for(int columnIndex=0;columnIndex<hObject.getSelectedMemberCount();columnIndex++)
                {
                    Datatype dt=hObject.getSelectedMemberTypes()[columnIndex];
                    //insertIntoCell((columnIndex+1), (index+1), hObject.getMemberNames()[columnIndex], xSpreadsheet, "F");
                }
                String[] rowData=(String[])o.elementAt(index);
                data[index]=new SVLVar(rowData);
//                for(int columnIndex=0;columnIndex<rowData.length;columnIndex++)
//                {
//                    //scores[columnIndex]=(double)rowData[columnIndex];
//                }
            }
            else
            {
                //double[] scores = new double[11];
                double[] rowData=(double[])o.elementAt(index);
                data[index]=new SVLVar(rowData);
//                for(int columnIndex=0;columnIndex<rowData.length;columnIndex++)
//                {
//                    scores[columnIndex]=(double)rowData[columnIndex];
//                }
            }
        }
            h5File.close();
    return new SVLVar(data);
}

private SVLVar retrieveResiduePWD(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
        String xPath="/DivCon/"+target+"/PWD Solvation";
        HObject ho=findHDF5Object(h5File, xPath);
        H5Group hObject=(H5Group)ho;
        if(hObject==null) throw new SVLJavaException("PWD does not exist for: '" + target + "'.");
    SVLVar[] data=new SVLVar[hObject.getNumberOfMembersInFile()];
        List pwdMembers=hObject.getMemberList();
        for(int index=0;index<hObject.getNumberOfMembersInFile();index++)
        {
            H5Group pwdGroup=(H5Group)pwdMembers.get(index);
            SVLVar[] pwdDataset=new SVLVar[pwdGroup.getNumberOfMembersInFile()+1];
            if(pwdDataset.length!=5) throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            List pwdRow=pwdGroup.getMemberList();
            pwdDataset[0]=new SVLVar(pwdGroup.getName());
            H5ScalarDS sequenceA=(H5ScalarDS)pwdRow.get(0);
            pwdDataset[1]=new SVLVar((int[])sequenceA.read());
            H5ScalarDS sequenceB=(H5ScalarDS)pwdRow.get(1);
            pwdDataset[2]=new SVLVar((int[])sequenceB.read());
            H5ScalarDS sequenceValues=(H5ScalarDS)pwdRow.get(2);
            pwdDataset[3]=new SVLVar((double[])sequenceValues.read());
            H5ScalarDS sequenceLabels=(H5ScalarDS)pwdRow.get(3);
            pwdDataset[4]=new SVLVar(sequenceLabels.readBytes());
            data[index]=new SVLVar(pwdDataset);
        }
            h5File.close();
    return new SVLVar(data);
}
    
private SVLVar retrieveAtomByAtomPWD(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
                    String xPath="/DivCon/"+target+"/PWD";
                    H5Group pwdGroup=(H5Group)findHDF5Object(h5File, xPath);
                    List<HObject> ligandList=pwdGroup.getMemberList();
                    H5Group pwdTargetLigandObject=(H5Group)ligandList.get(0);
                    H5ScalarDS pwdObject=(H5ScalarDS)pwdTargetLigandObject.getMemberList().get(0);
                    double[] pwdValues= (double[])pwdObject.read();  
                    FileWriter pwdOutputFile=new FileWriter("./"+target+"-"+pwdObject.getName()+".pwd");
//                    BufferedWriter bufferedPWDWriter=new BufferedWriter(pwdOutputFile);
                    xPath="/DivCon/"+target+"/"+target;
                    H5CompoundDS hTargetCollectionObject=(H5CompoundDS)findHDF5Object(h5File, xPath);
                    hTargetCollectionObject.selectMember(0);
                    String[] targetAtomSymbols=(String[])((Vector)hTargetCollectionObject.read()).elementAt(0);
                    hTargetCollectionObject.selectMember(1);
                    String[] targetAtomNames=(String[])((Vector)hTargetCollectionObject.read()).elementAt(1);
                    hTargetCollectionObject.selectMember(2);
                    String[] targetResidueNames=(String[])((Vector)hTargetCollectionObject.read()).elementAt(2);
                    hTargetCollectionObject.selectMember(3);
                    int[] targetSequence=(int[])((Vector)hTargetCollectionObject.read()).elementAt(3);
                    hTargetCollectionObject.selectMember(4);
                    double[] x=(double[])((Vector)hTargetCollectionObject.read()).elementAt(4);
                    hTargetCollectionObject.selectMember(5);
                    double[] y=(double[])((Vector)hTargetCollectionObject.read()).elementAt(5);
                    hTargetCollectionObject.selectMember(6);
                    double[] z=(double[])((Vector)hTargetCollectionObject.read()).elementAt(6);
                    hTargetCollectionObject.selectMember(7);
                    double[] epsilon=(double[])((Vector)hTargetCollectionObject.read()).elementAt(7);
                    hTargetCollectionObject.selectMember(8);
                    double[] welllDepth=(double[])((Vector)hTargetCollectionObject.read()).elementAt(8);
                    hTargetCollectionObject.selectMember(9);
                    int[] formalCharge=(int[])((Vector)hTargetCollectionObject.read()).elementAt(9);
                    xPath="/DivCon/"+target+"/Trajectories/Geometries/Final";
                    H5ScalarDS trajectoryObject=(H5ScalarDS)findHDF5Object(h5File, xPath);
                    SVLJava.print(" trajectoryObject? "+trajectoryObject);
                    if(trajectoryObject!=null)
                    {
                    double[] buf= (double[])trajectoryObject.read();  
                    long[] dims=trajectoryObject.getDims();
                    for(int i=0;i<dims[0];i++)
                    {
                    x[i]=buf[i*3];
                    y[i]=buf[i*3+1];
                    z[i]=buf[i*3+2];
                    }
                    }
                    ListIterator li=hTargetCollectionObject.getMetadata().listIterator();
                    while(li.hasNext())
                    {
                        Object obj=li.next();
                        if(obj instanceof ncsa.hdf.object.Attribute && ((ncsa.hdf.object.Attribute)obj).getName().compareTo("Total Charge")==0)
                        {
//                            ChargesType chargeType=new ChargesType();
//                            TotalChargeType totalChargeType=new TotalChargeType();
//                            chargeType.setTotalCharge(totalChargeType);
//                            int[] tc=(int[])((ncsa.hdf.object.Attribute)obj).getValue();
                        }
                    }
                    xPath="/DivCon/"+pwdObject.getName();
                    H5CompoundDS hLigandCollectionObject=(H5CompoundDS)findHDF5Object(h5File, xPath);
                    hLigandCollectionObject.selectMember(0);
                    String[] ligandAtomSymbols=(String[])((Vector)hLigandCollectionObject.read()).elementAt(0);
                    hLigandCollectionObject.selectMember(1);
                    String[] ligandAtomNames=(String[])((Vector)hLigandCollectionObject.read()).elementAt(1);
                    hLigandCollectionObject.selectMember(2);
                    String[] ligandResidueNames=(String[])((Vector)hLigandCollectionObject.read()).elementAt(2);
                    hLigandCollectionObject.selectMember(3);
                    int[] ligandSequence=(int[])((Vector)hLigandCollectionObject.read()).elementAt(3);
                    hLigandCollectionObject.selectMember(4);
                    double[] xLigand=(double[])((Vector)hLigandCollectionObject.read()).elementAt(4);
                    hLigandCollectionObject.selectMember(5);
                    double[] yLigand=(double[])((Vector)hLigandCollectionObject.read()).elementAt(5);
                    hLigandCollectionObject.selectMember(6);
                    double[] zLigand=(double[])((Vector)hLigandCollectionObject.read()).elementAt(6);
//                    hLigandCollectionObject.selectMember(7);
//                    double[] epsilon=(double[])((Vector)hLigandCollectionObject.read()).elementAt(7);
//                    hLigandCollectionObject.selectMember(8);
//                    double[] welllDepth=(double[])((Vector)hLigandCollectionObject.read()).elementAt(8);
//                    hLigandCollectionObject.selectMember(9);
//                    int[] formalCharge=(int[])((Vector)hLigandCollectionObject.read()).elementAt(9);
                    SVLJava.print(targetAtomSymbols.length+" lengths "+ligandAtomSymbols.length);
    SVLVar[] data=new SVLVar[10];
    ArrayList<Integer> targetIndexList=new  ArrayList<Integer>();
    ArrayList<Integer> ligandIndexList=new  ArrayList<Integer>();
    ArrayList<Double> EabList=new  ArrayList<Double>();
    ArrayList<Double> EabpList=new  ArrayList<Double>();
    ArrayList<Double> EabcList=new  ArrayList<Double>();
    ArrayList<Double> LennardJonesList=new  ArrayList<Double>();
    ArrayList<Double> dispersionList=new  ArrayList<Double>();
    ArrayList<Double> repulsionList=new  ArrayList<Double>();
    ArrayList<Double> electrostaticList=new  ArrayList<Double>();
    ArrayList<Double> distanceList=new  ArrayList<Double>();
                        //bufferedPWDWriter.append("Target Index,Target Atom,Target Residue,Target Sequence,Ligand Index,Ligand Atom,Ligand Residue,Ligand Sequence,Distance,E_AB,E_ABp,E_ABc\n");
                    for(int atomIndex=0;atomIndex<targetAtomSymbols.length;atomIndex++)
                    {
                    for(int ligandAtomIndex=0;ligandAtomIndex<ligandAtomSymbols.length;ligandAtomIndex++)
                    {
                        int iw=indexW(atomIndex+1, targetAtomSymbols.length+ligandAtomIndex+1);
                        double distance=Math.sqrt(Math.pow(x[atomIndex]-xLigand[ligandAtomIndex],2)+Math.pow(y[atomIndex]-yLigand[ligandAtomIndex],2)+Math.pow(z[atomIndex]-zLigand[ligandAtomIndex],2));
                        if(distance<=4.0 || pwdValues[7*iw]!=0.0)
                        {
                            targetIndexList.add(atomIndex);
                            ligandIndexList.add(ligandAtomIndex);
                            EabList.add(pwdValues[7*iw]);
                            EabpList.add(pwdValues[7*iw+1]);
                            EabcList.add(pwdValues[7*iw+2]);
                            LennardJonesList.add(pwdValues[7*iw+3]);
                            dispersionList.add(pwdValues[7*iw+4]);
                            repulsionList.add(pwdValues[7*iw+5]);
                            electrostaticList.add(pwdValues[7*iw+6]);
                            distanceList.add(distance);
//                            //bufferedPWDWriter.append(atomIndex+","+targetAtomNames[atomIndex]+","+targetResidueNames[atomIndex]+","+targetSequence[atomIndex]+","+ligandAtomIndex+","+ligandAtomNames[ligandAtomIndex]+","+ligandResidueNames[ligandAtomIndex]+","+ligandSequence[ligandAtomIndex]+","+distance+","+pwdValues[3*iw]+","+pwdValues[3*iw+1]+","+pwdValues[3*iw+2]+"\n");
                        }
                    }
                    }
                    int[] targetIndex=new int[targetIndexList.size()];
                    int[] ligandIndex=new int[targetIndexList.size()];
                    double[] Eab=new double[targetIndexList.size()];
                    double[] Eabp=new double[targetIndexList.size()];
                    double[] Eabc=new double[targetIndexList.size()];
                    double[] LennardJones=new double[targetIndexList.size()];
                    double[] dispersion=new double[targetIndexList.size()];
                    double[] repulsion=new double[targetIndexList.size()];
                    double[] electrostatic=new double[targetIndexList.size()];
                    double[] distance=new double[targetIndexList.size()];
                    for(int index=0;index<targetIndex.length;index++)
                    {
                    targetIndex[index]=targetIndexList.get(index);
                    ligandIndex[index]=ligandIndexList.get(index);
                    Eab[index]=EabList.get(index);
                    Eabp[index]=EabpList.get(index);
                    Eabc[index]=EabcList.get(index);
                    LennardJones[index]=LennardJonesList.get(index);
                    dispersion[index]=dispersionList.get(index);
                    repulsion[index]=repulsionList.get(index);
                    electrostatic[index]=electrostaticList.get(index);
                    distance[index]=distanceList.get(index);
                    }
    data[0]=new SVLVar(targetIndex);
    data[1]=new SVLVar(ligandIndex);
    data[2]=new SVLVar(Eab);
    data[3]=new SVLVar(Eabp);
    data[4]=new SVLVar(Eabc);
    data[5]=new SVLVar(LennardJones);
    data[6]=new SVLVar(dispersion);
    data[7]=new SVLVar(repulsion);
    data[8]=new SVLVar(electrostatic);
    data[9]=new SVLVar(distance);
    
                    
                    //bufferedPWDWriter.close();
//        List pwdMembers=hObject.getMemberList();
//        for(int index=0;index<hObject.getNumberOfMembersInFile();index++)
//        {
//            H5Group pwdGroup=(H5Group)pwdMembers.get(index);
//            SVLVar[] pwdDataset=new SVLVar[pwdGroup.getNumberOfMembersInFile()+1];
//            if(pwdDataset.length!=5) throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
//            List pwdRow=pwdGroup.getMemberList();
//            pwdDataset[0]=new SVLVar(pwdGroup.getName());
//            H5ScalarDS sequenceA=(H5ScalarDS)pwdRow.get(0);
//            pwdDataset[1]=new SVLVar((int[])sequenceA.read());
//            H5ScalarDS sequenceB=(H5ScalarDS)pwdRow.get(1);
//            pwdDataset[2]=new SVLVar((int[])sequenceB.read());
//            H5ScalarDS sequenceValues=(H5ScalarDS)pwdRow.get(2);
//            pwdDataset[3]=new SVLVar((double[])sequenceValues.read());
//            H5ScalarDS sequenceLabels=(H5ScalarDS)pwdRow.get(3);
//            pwdDataset[4]=new SVLVar(sequenceLabels.readBytes());
//            data[index]=new SVLVar(pwdDataset);
//        }
            h5File.close();
    return new SVLVar(data);
}
    
  protected HObject findHDF5Object(FileFormat file, String path)
  {
      if (file == null || path == null)
          return null;

      if (!path.endsWith("/"))
          path = path+"/";

      DefaultMutableTreeNode theRoot =
                        (DefaultMutableTreeNode)file.getRootNode();

      if (theRoot == null)
          return null;
      else if (path.equals("/"))
          return (HObject)theRoot.getUserObject();

      Enumeration local_enum =
            ((DefaultMutableTreeNode)theRoot).breadthFirstEnumeration();
      DefaultMutableTreeNode theNode = null;
      HObject theObj = null;
      boolean hit=false;
      while(local_enum.hasMoreElements())
      {
          theNode = (DefaultMutableTreeNode)local_enum.nextElement();
          theObj = (HObject)theNode.getUserObject();
          String fullPath = theObj.getFullName()+"/";
          if (path.equals(fullPath))
          {
              hit=true;
              break;
          }
      }
      if(!hit) return null;
      return theObj;
 }
 
    protected int indexW(int i, int j)
    {
        int indw,ki,kj;
        ki=Math.max(i,j);
        kj=Math.min(i,j);
        indw=(ki*(ki-1))/2+kj;
        return indw;
    }
    
    
}
