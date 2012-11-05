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
import com.quantumbioinc.datacorrespondent.Correspondent;

        
/**
 *
 * @author roger
 */
public class HDF5Correspondent extends Correspondent implements SVLJavaDriver {

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
	    if(cmd.equals("listModels" ))
            {
                res = listModels(res);
            }
            else if(cmd.equals("storeModel" ))
            {
                res = storeModel(res);
            }
            else if(cmd.equals("retrieveQMScore" ))
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
	    else if (cmd.equals("retrieveAtomByAtomMRM"))
            {
                res = retrieveAtomByAtomMRM(res);
            }
	    else if (cmd.equals("retrieveNMRAverages"))
            {
                res = retrieveNMRAverages(res);
            }
	    else {
		throw new SVLJavaException("here Unknown command: '" + cmd + "'.");
	    }
	    return res;
	}
	catch (Exception ex) {
	    StringWriter sw = new StringWriter();
	    ex.printStackTrace(new PrintWriter(sw));
	    throw new SVLJavaException("Exception:\n"+sw.toString());
	}
    }

    private SVLVar listModels(SVLVar var) throws SVLJavaException, IOException, Exception
    {
        String filename = var.peek(1).getTokn(1);
        H5File h5File=new H5File(filename, H5File.READ);
        h5File.open();
        ArrayList<String> targetList=findTargets(h5File);
        SVLVar[] data=new SVLVar[targetList.size()];
        for(int index=0;index<targetList.size();index++)
        {
            SVLVar[] individuals=new SVLVar[2];
            individuals[0]=new SVLVar(targetList.get(index), true);
            //individuals[0].
            //individuals[1]=new SVLVar("Target");
            String xPath="/DivCon/"+targetList.get(index)+"/QM Score";
            HObject ho=findHDF5Object(h5File, xPath);
            H5CompoundDS hObject=(H5CompoundDS)ho;
            if(hObject!=null){
                Vector o=(Vector)hObject.getData();
                individuals[1]=new SVLVar((String[])o.elementAt(0));
            }
            else
            {
                individuals[1]=new SVLVar();
            }
            data[index]=new SVLVar(new String[]{"target", "ligand"}, individuals);
        }
        h5File.close();
        return new SVLVar(data);
    }

    private SVLVar storeModel(SVLVar var) throws SVLJavaException, IOException, Exception
    {
        String filename = var.peek(1).getTokn(1);
        String target = var.peek(1).getTokn(2);
        H5File h5File=new H5File(filename, H5File.READ);
        h5File.open();
                String xPath="/DivCon";
                H5Group divconGroup=(H5Group)findHDF5Object(h5File, xPath);
                H5Group targetGroup=new H5Group(h5File, target, xPath, divconGroup);
                divconGroup.addToMemberList(targetGroup);
        h5File.close();
        return new SVLVar();
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
        if(hObject==null) return new SVLVar();//throw new SVLJavaException("PWD does not exist for: '" + target + "'.");
    SVLVar[] data=new SVLVar[hObject.getNumberOfMembersInFile()];
        List pwdMembers=hObject.getMemberList();
        for(int index=0;index<hObject.getNumberOfMembersInFile();index++)
        {
            H5Group pwdGroup=(H5Group)pwdMembers.get(index);
            SVLVar[] pwdDataset=new SVLVar[pwdGroup.getNumberOfMembersInFile()+1];
            if(pwdDataset.length!=5) throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            List pwdRow=pwdGroup.getMemberList();
            pwdDataset[0]=new SVLVar(pwdGroup.getName());
            for(int memberIndex=0;memberIndex<4;memberIndex++)
            {
                H5ScalarDS member=(H5ScalarDS)pwdRow.get(memberIndex);
                if(member.getName().compareTo("Sequence A")==0)
                {
                    pwdDataset[1]=new SVLVar((int[])member.read());
                }
                else if(member.getName().compareTo("Sequence B")==0)
                {
                    pwdDataset[2]=new SVLVar((int[])member.read());
                }
                else if(member.getName().compareTo(target)==0)
                {
                    pwdDataset[4]=new SVLVar(member.readBytes());
                }
                else
                {
                    pwdDataset[3]=new SVLVar((double[])member.read());
                }
            }
//            pwdDataset[1]=new SVLVar((int[])sequenceA.read());
//            H5ScalarDS sequenceB=(H5ScalarDS)pwdRow.get(1);
//            pwdDataset[2]=new SVLVar((int[])sequenceB.read());
//            H5ScalarDS sequenceValues=(H5ScalarDS)pwdRow.get(2);
//            pwdDataset[3]=new SVLVar((double[])sequenceValues.read());
//            H5ScalarDS sequenceLabels=(H5ScalarDS)pwdRow.get(3);
//            pwdDataset[4]=new SVLVar(sequenceLabels.readBytes());
            data[index]=new SVLVar(new String[]{"name", "indexA", "indexB", "values", "labels"}, pwdDataset);
        }
            h5File.close();
    return new SVLVar(data);
}
    
private SVLVar retrieveAtomByAtomPWD(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    String ligand="";
    if(var.peek(1).length()>2)
    {
        ligand = var.peek(1).getTokn(3);
        //if(true) return new SVLVar("retrieveAtomByAtomPWD "+ligand);
    }
//                    FileWriter pwdOutputFile=new FileWriter("./check.pwd");
//                    BufferedWriter bufferedPWDWriter=new BufferedWriter(pwdOutputFile);
//                        bufferedPWDWriter.append(target+" here\n");
//                        bufferedPWDWriter.flush();
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
                    String xPath="/DivCon/"+target+"/PWD";
                    H5Group pwdGroup=(H5Group)findHDF5Object(h5File, xPath);
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
    if(ligand.length()>0)
    {
            List<HObject> ligandList=pwdGroup.getMemberList();
            H5Group pwdTargetLigandObject=null;
            boolean hit=false;
            for(int index=0;index<pwdGroup.getNumberOfMembersInFile();index++)
            {
                pwdTargetLigandObject=(H5Group)ligandList.get(index);
                if(pwdTargetLigandObject.getName().endsWith("-"+ligand))
                {
                    hit=true;
                    break;
                }
            }
            if(!hit)
            {
                h5File.close();
            return new SVLVar();
            }
            SVLVar[] pwdDataset=new SVLVar[pwdTargetLigandObject.getNumberOfMembersInFile()+8];
            if(pwdDataset.length!=11) throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            List pwdRow=pwdTargetLigandObject.getMemberList();
            pwdDataset[0]=new SVLVar(pwdTargetLigandObject.getName());
            int[] targetIndex=null;
            int[] ligandIndex=null;
            double[] pwdValues=null;
            H5ScalarDS pwdObject=null;
            for(int memberIndex=0;memberIndex<3;memberIndex++)
            {
                H5ScalarDS member=(H5ScalarDS)pwdRow.get(memberIndex);
                if(member.getName().compareTo("Index A")==0)
                {
                    targetIndex=(int[])member.read();
                    pwdDataset[1]=new SVLVar((int[])member.read());
                }
                else if(member.getName().compareTo("Index B")==0)
                {
                    ligandIndex=(int[])member.read();
                    pwdDataset[2]=new SVLVar((int[])member.read());
                }
                else
                {
                    pwdObject=(H5ScalarDS)member;
                    pwdValues=(double[])member.read();
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
                    ArrayList<Double> EabList=new  ArrayList<Double>();
                    ArrayList<Double> EabpList=new  ArrayList<Double>();
                    ArrayList<Double> EabcList=new  ArrayList<Double>();
                    ArrayList<Double> LennardJonesList=new  ArrayList<Double>();
                    ArrayList<Double> dispersionList=new  ArrayList<Double>();
                    ArrayList<Double> repulsionList=new  ArrayList<Double>();
                    ArrayList<Double> electrostaticList=new  ArrayList<Double>();
                    ArrayList<Double> distanceList=new  ArrayList<Double>();
                    
//                        bufferedPWDWriter.append(targetIndex.length+" "+ligandIndex.length+"\n");
                    for(int atomIndex=0;atomIndex<targetIndex.length;atomIndex++)
                    {
                        int ligandOffsetIndex=ligandIndex[atomIndex];
//                        bufferedPWDWriter.append(targetIndex[atomIndex]+" "+ligandOffsetIndex+"\n");
                        ligandOffsetIndex-=x.length;
                        double distance=Math.sqrt(Math.pow(x[targetIndex[atomIndex]]-xLigand[ligandOffsetIndex],2)+Math.pow(y[targetIndex[atomIndex]]-yLigand[ligandOffsetIndex],2)+Math.pow(z[targetIndex[atomIndex]]-zLigand[ligandOffsetIndex],2));
                        //if(distance<=4.0 || pwdValues[7*iw]!=0.0)
                        EabList.add(pwdValues[7*atomIndex]);
                        EabpList.add(pwdValues[7*atomIndex+1]);
                        EabcList.add(pwdValues[7*atomIndex+2]);
                        LennardJonesList.add(pwdValues[7*atomIndex+3]);
                        dispersionList.add(pwdValues[7*atomIndex+4]);
                        repulsionList.add(pwdValues[7*atomIndex+5]);
                        electrostaticList.add(pwdValues[7*atomIndex+6]);
                        distanceList.add(distance);
                    }
                    double[] Eab=new double[pwdDataset[1].length()];
                    double[] Eabp=new double[pwdDataset[1].length()];
                    double[] Eabc=new double[pwdDataset[1].length()];
                    double[] LennardJones=new double[pwdDataset[1].length()];
                    double[] dispersion=new double[pwdDataset[1].length()];
                    double[] repulsion=new double[pwdDataset[1].length()];
                    double[] electrostatic=new double[pwdDataset[1].length()];
                    double[] distance=new double[pwdDataset[1].length()];
                    for(int vIndex=0;vIndex<pwdDataset[1].length();vIndex++)
                    {
                    Eab[vIndex]=EabList.get(vIndex);
                    Eabp[vIndex]=EabpList.get(vIndex);
                    Eabc[vIndex]=EabcList.get(vIndex);
                    LennardJones[vIndex]=LennardJonesList.get(vIndex);
                    dispersion[vIndex]=dispersionList.get(vIndex);
                    repulsion[vIndex]=repulsionList.get(vIndex);
                    electrostatic[vIndex]=electrostaticList.get(vIndex);
                    distance[vIndex]=distanceList.get(vIndex);
                    }
                    pwdDataset[3]=new SVLVar(Eab);
                    pwdDataset[4]=new SVLVar(Eabp);
                    pwdDataset[5]=new SVLVar(Eabc);
                    pwdDataset[6]=new SVLVar(LennardJones);
                    pwdDataset[7]=new SVLVar(dispersion);
                    pwdDataset[8]=new SVLVar(repulsion);
                    pwdDataset[9]=new SVLVar(electrostatic);
                    pwdDataset[10]=new SVLVar(distance);
            SVLVar data=new SVLVar(new String[]{"name", "indexA", "indexB", "Eab", "Eabp", "Eabc", "LennardJones", "dispersion", "repulsion", "electrostatic", "distance"}, pwdDataset);
        h5File.close();
    return new SVLVar(data);
    }
    else
    {
    SVLVar[] data=new SVLVar[pwdGroup.getNumberOfMembersInFile()];
    
        for(int index=0;index<pwdGroup.getNumberOfMembersInFile();index++)
        {
                    List<HObject> ligandList=pwdGroup.getMemberList();
                    H5Group pwdTargetLigandObject=(H5Group)ligandList.get(index);
            SVLVar[] pwdDataset=new SVLVar[pwdTargetLigandObject.getNumberOfMembersInFile()+8];
            if(pwdDataset.length!=11) throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            List pwdRow=pwdTargetLigandObject.getMemberList();
            pwdDataset[0]=new SVLVar(pwdTargetLigandObject.getName());
            int[] targetIndex=null;
            int[] ligandIndex=null;
            double[] pwdValues=null;
            H5ScalarDS pwdObject=null;
            for(int memberIndex=0;memberIndex<3;memberIndex++)
            {
                H5ScalarDS member=(H5ScalarDS)pwdRow.get(memberIndex);
                if(member.getName().compareTo("Index A")==0)
                {
                    targetIndex=(int[])member.read();
                    pwdDataset[1]=new SVLVar((int[])member.read());
                }
                else if(member.getName().compareTo("Index B")==0)
                {
                    ligandIndex=(int[])member.read();
                    pwdDataset[2]=new SVLVar((int[])member.read());
                }
                else
                {
                    pwdObject=(H5ScalarDS)member;
                    pwdValues=(double[])member.read();
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
                    ArrayList<Double> EabList=new  ArrayList<Double>();
                    ArrayList<Double> EabpList=new  ArrayList<Double>();
                    ArrayList<Double> EabcList=new  ArrayList<Double>();
                    ArrayList<Double> LennardJonesList=new  ArrayList<Double>();
                    ArrayList<Double> dispersionList=new  ArrayList<Double>();
                    ArrayList<Double> repulsionList=new  ArrayList<Double>();
                    ArrayList<Double> electrostaticList=new  ArrayList<Double>();
                    ArrayList<Double> distanceList=new  ArrayList<Double>();
                    
//                        bufferedPWDWriter.append(targetIndex.length+" "+ligandIndex.length+"\n");
                    for(int atomIndex=0;atomIndex<targetIndex.length;atomIndex++)
                    {
                        int ligandOffsetIndex=ligandIndex[atomIndex];
//                        bufferedPWDWriter.append(targetIndex[atomIndex]+" "+ligandOffsetIndex+"\n");
                        ligandOffsetIndex-=x.length;
                        double distance=Math.sqrt(Math.pow(x[targetIndex[atomIndex]]-xLigand[ligandOffsetIndex],2)+Math.pow(y[targetIndex[atomIndex]]-yLigand[ligandOffsetIndex],2)+Math.pow(z[targetIndex[atomIndex]]-zLigand[ligandOffsetIndex],2));
                        //if(distance<=4.0 || pwdValues[7*iw]!=0.0)
                        EabList.add(pwdValues[7*atomIndex]);
                        EabpList.add(pwdValues[7*atomIndex+1]);
                        EabcList.add(pwdValues[7*atomIndex+2]);
                        LennardJonesList.add(pwdValues[7*atomIndex+3]);
                        dispersionList.add(pwdValues[7*atomIndex+4]);
                        repulsionList.add(pwdValues[7*atomIndex+5]);
                        electrostaticList.add(pwdValues[7*atomIndex+6]);
                        distanceList.add(distance);
                    }
                    double[] Eab=new double[pwdDataset[1].length()];
                    double[] Eabp=new double[pwdDataset[1].length()];
                    double[] Eabc=new double[pwdDataset[1].length()];
                    double[] LennardJones=new double[pwdDataset[1].length()];
                    double[] dispersion=new double[pwdDataset[1].length()];
                    double[] repulsion=new double[pwdDataset[1].length()];
                    double[] electrostatic=new double[pwdDataset[1].length()];
                    double[] distance=new double[pwdDataset[1].length()];
                    for(int vIndex=0;vIndex<pwdDataset[1].length();vIndex++)
                    {
                    Eab[vIndex]=EabList.get(vIndex);
                    Eabp[vIndex]=EabpList.get(vIndex);
                    Eabc[vIndex]=EabcList.get(vIndex);
                    LennardJones[vIndex]=LennardJonesList.get(vIndex);
                    dispersion[vIndex]=dispersionList.get(vIndex);
                    repulsion[vIndex]=repulsionList.get(vIndex);
                    electrostatic[vIndex]=electrostaticList.get(vIndex);
                    distance[vIndex]=distanceList.get(vIndex);
                    }
                    pwdDataset[3]=new SVLVar(Eab);
                    pwdDataset[4]=new SVLVar(Eabp);
                    pwdDataset[5]=new SVLVar(Eabc);
                    pwdDataset[6]=new SVLVar(LennardJones);
                    pwdDataset[7]=new SVLVar(dispersion);
                    pwdDataset[8]=new SVLVar(repulsion);
                    pwdDataset[9]=new SVLVar(electrostatic);
                    pwdDataset[10]=new SVLVar(distance);
            data[index]=new SVLVar(new String[]{"name", "indexA", "indexB", "Eab", "Eabp", "Eabc", "LennardJones", "dispersion", "repulsion", "electrostatic", "distance"}, pwdDataset);
        }
//             bufferedPWDWriter.close();       
        h5File.close();
    return new SVLVar(data);
    }
}
    
private SVLVar retrieveAtomByAtomMRM(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    String ligand="";
    if(var.peek(1).length()>2)
    {
        ligand = var.peek(1).getTokn(3);
        //if(true) return new SVLVar("retrieveAtomByAtomMRM "+ligand);
    }
//                    FileWriter pwdOutputFile=new FileWriter("./check.pwd");
//                    BufferedWriter bufferedPWDWriter=new BufferedWriter(pwdOutputFile);
//                        bufferedPWDWriter.append(target+" here\n");
//                        bufferedPWDWriter.flush();
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
                    String xPath="/DivCon/"+target+"/MRM";
                    H5Group pwdGroup=(H5Group)findHDF5Object(h5File, xPath);
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
    if(ligand.length()>0)
    {
            List<HObject> ligandList=pwdGroup.getMemberList();
            H5Group pwdTargetLigandObject=null;
            boolean hit=false;
            for(int index=0;index<pwdGroup.getNumberOfMembersInFile();index++)
            {
                pwdTargetLigandObject=(H5Group)ligandList.get(index);
                if(pwdTargetLigandObject.getName().endsWith("-"+ligand))
                {
                    hit=true;
                    break;
                }
            }
            if(!hit)
            {
                h5File.close();
            return new SVLVar();
            }
            SVLVar[] pwdDataset=new SVLVar[pwdTargetLigandObject.getNumberOfMembersInFile()+4];
            if(pwdDataset.length!=7) throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            List pwdRow=pwdTargetLigandObject.getMemberList();
            pwdDataset[0]=new SVLVar(pwdTargetLigandObject.getName());
            int[] targetIndex=null;
            int[] ligandIndex=null;
            double[] pwdValues=null;
            H5ScalarDS pwdObject=null;
            for(int memberIndex=0;memberIndex<3;memberIndex++)
            {
                H5ScalarDS member=(H5ScalarDS)pwdRow.get(memberIndex);
                if(member.getName().compareTo("Index A")==0)
                {
                    targetIndex=(int[])member.read();
                    pwdDataset[1]=new SVLVar((int[])member.read());
                }
                else if(member.getName().compareTo("Index B")==0)
                {
                    ligandIndex=(int[])member.read();
                    pwdDataset[2]=new SVLVar((int[])member.read());
                }
                else
                {
                    pwdObject=(H5ScalarDS)member;
                    pwdValues=(double[])member.read();
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
                    ArrayList<Double> EabList=new  ArrayList<Double>();
                    ArrayList<Double> EabpList=new  ArrayList<Double>();
                    ArrayList<Double> EabcList=new  ArrayList<Double>();
                    ArrayList<Double> LennardJonesList=new  ArrayList<Double>();
                    ArrayList<Double> dispersionList=new  ArrayList<Double>();
                    ArrayList<Double> repulsionList=new  ArrayList<Double>();
                    ArrayList<Double> electrostaticList=new  ArrayList<Double>();
                    ArrayList<Double> distanceList=new  ArrayList<Double>();
                    
//                        bufferedPWDWriter.append(targetIndex.length+" "+ligandIndex.length+"\n");
                    for(int atomIndex=0;atomIndex<targetIndex.length;atomIndex++)
                    {
                        int ligandOffsetIndex=ligandIndex[atomIndex];
//                        bufferedPWDWriter.append(targetIndex[atomIndex]+" "+ligandOffsetIndex+"\n");
                        double distance=Math.sqrt(Math.pow(x[targetIndex[atomIndex]]-xLigand[ligandOffsetIndex],2)+Math.pow(y[targetIndex[atomIndex]]-yLigand[ligandOffsetIndex],2)+Math.pow(z[targetIndex[atomIndex]]-zLigand[ligandOffsetIndex],2));
                        //if(distance<=4.0 || pwdValues[7*iw]!=0.0)
                        EabList.add(pwdValues[3*atomIndex]);
                        EabpList.add(pwdValues[3*atomIndex+1]);
                        EabcList.add(pwdValues[3*atomIndex+2]);
                        distanceList.add(distance);
                    }
                    double[] Eab=new double[pwdDataset[1].length()];
                    double[] Eabp=new double[pwdDataset[1].length()];
                    double[] Eabc=new double[pwdDataset[1].length()];
                    double[] distance=new double[pwdDataset[1].length()];
                    for(int vIndex=0;vIndex<pwdDataset[1].length();vIndex++)
                    {
                    Eab[vIndex]=EabList.get(vIndex);
                    Eabp[vIndex]=EabpList.get(vIndex);
                    Eabc[vIndex]=EabcList.get(vIndex);
                    distance[vIndex]=distanceList.get(vIndex);
                    }
                    pwdDataset[3]=new SVLVar(Eab);
                    pwdDataset[4]=new SVLVar(Eabp);
                    pwdDataset[5]=new SVLVar(Eabc);
                    pwdDataset[6]=new SVLVar(distance);
        SVLVar data=new SVLVar(new String[]{"name", "indexA", "indexB", "mrmSteric", "mrmElectrostatic", "mrmPhenomenologicalSolvation", "distance"}, pwdDataset);
        h5File.close();
    return new SVLVar(data);
    }
    else
    {
        SVLVar[] data=new SVLVar[pwdGroup.getNumberOfMembersInFile()];
    
        for(int index=0;index<pwdGroup.getNumberOfMembersInFile();index++)
        {
                    List<HObject> ligandList=pwdGroup.getMemberList();
                    H5Group pwdTargetLigandObject=(H5Group)ligandList.get(index);
            SVLVar[] pwdDataset=new SVLVar[pwdTargetLigandObject.getNumberOfMembersInFile()+4];
            if(pwdDataset.length!=7) throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            List pwdRow=pwdTargetLigandObject.getMemberList();
            pwdDataset[0]=new SVLVar(pwdTargetLigandObject.getName());
            int[] targetIndex=null;
            int[] ligandIndex=null;
            double[] pwdValues=null;
            H5ScalarDS pwdObject=null;
            for(int memberIndex=0;memberIndex<3;memberIndex++)
            {
                H5ScalarDS member=(H5ScalarDS)pwdRow.get(memberIndex);
                if(member.getName().compareTo("Index A")==0)
                {
                    targetIndex=(int[])member.read();
                    pwdDataset[1]=new SVLVar((int[])member.read());
                }
                else if(member.getName().compareTo("Index B")==0)
                {
                    ligandIndex=(int[])member.read();
                    pwdDataset[2]=new SVLVar((int[])member.read());
                }
                else
                {
                    pwdObject=(H5ScalarDS)member;
                    pwdValues=(double[])member.read();
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
                    ArrayList<Double> EabList=new  ArrayList<Double>();
                    ArrayList<Double> EabpList=new  ArrayList<Double>();
                    ArrayList<Double> EabcList=new  ArrayList<Double>();
                    ArrayList<Double> LennardJonesList=new  ArrayList<Double>();
                    ArrayList<Double> dispersionList=new  ArrayList<Double>();
                    ArrayList<Double> repulsionList=new  ArrayList<Double>();
                    ArrayList<Double> electrostaticList=new  ArrayList<Double>();
                    ArrayList<Double> distanceList=new  ArrayList<Double>();
                    
//                        bufferedPWDWriter.append(targetIndex.length+" "+ligandIndex.length+"\n");
                    for(int atomIndex=0;atomIndex<targetIndex.length;atomIndex++)
                    {
                        int ligandOffsetIndex=ligandIndex[atomIndex];
//                        bufferedPWDWriter.append(targetIndex[atomIndex]+" "+ligandOffsetIndex+"\n");
                        double distance=Math.sqrt(Math.pow(x[targetIndex[atomIndex]]-xLigand[ligandOffsetIndex],2)+Math.pow(y[targetIndex[atomIndex]]-yLigand[ligandOffsetIndex],2)+Math.pow(z[targetIndex[atomIndex]]-zLigand[ligandOffsetIndex],2));
                        //if(distance<=4.0 || pwdValues[7*iw]!=0.0)
                        EabList.add(pwdValues[3*atomIndex]);
                        EabpList.add(pwdValues[3*atomIndex+1]);
                        EabcList.add(pwdValues[3*atomIndex+2]);
                        distanceList.add(distance);
                    }
                    double[] Eab=new double[pwdDataset[1].length()];
                    double[] Eabp=new double[pwdDataset[1].length()];
                    double[] Eabc=new double[pwdDataset[1].length()];
                    double[] distance=new double[pwdDataset[1].length()];
                    for(int vIndex=0;vIndex<pwdDataset[1].length();vIndex++)
                    {
                    Eab[vIndex]=EabList.get(vIndex);
                    Eabp[vIndex]=EabpList.get(vIndex);
                    Eabc[vIndex]=EabcList.get(vIndex);
                    distance[vIndex]=distanceList.get(vIndex);
                    }
                    pwdDataset[3]=new SVLVar(Eab);
                    pwdDataset[4]=new SVLVar(Eabp);
                    pwdDataset[5]=new SVLVar(Eabc);
                    pwdDataset[6]=new SVLVar(distance);
            data[index]=new SVLVar(new String[]{"name", "indexA", "indexB", "mrmSteric", "mrmElectrostatic", "mrmPhenomenologicalSolvation", "distance"}, pwdDataset);
        }
//             bufferedPWDWriter.close();       
        h5File.close();
    return new SVLVar(data);
    }
}
    
private SVLVar retrieveNMRAverages(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
                    String xPath="/DivCon/"+target+"/NMR Matrices";
                    H5Group nmrGroup=(H5Group)findHDF5Object(h5File, xPath);
            SVLVar[] averages=new SVLVar[3];
            if(nmrGroup!=null)
            {
                List nmrRow=nmrGroup.getMemberList();
                for(int memberIndex=0;memberIndex<5;memberIndex++)
                {
                    H5ScalarDS member=(H5ScalarDS)nmrRow.get(memberIndex);
                    if(member.getName().compareTo("Selected Indices")==0)
                    {
                        averages[0]=new SVLVar((int[])member.read());
                    }
                    else if(member.getName().compareTo("Average")==0)
                    {
                        averages[1]=new SVLVar((double[])member.read());
                    }
                    else if(member.getName().compareTo("Anisotropy")==0)
                    {
                        averages[2]=new SVLVar((double[])member.read());
                    }
                }
            }
            else
            {
                averages[0]=new SVLVar();
                averages[1]=new SVLVar();
                averages[2]=new SVLVar();
            }
    SVLVar data=new SVLVar(new String[]{"index", "average", "anisotropy"}, averages);
            h5File.close();
    return new SVLVar(data);
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
