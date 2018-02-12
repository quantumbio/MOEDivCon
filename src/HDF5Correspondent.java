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
import hdf.object.Attribute;
import hdf.object.Datatype;
import hdf.object.FileFormat;
import hdf.object.HObject;
import hdf.object.h5.H5CompoundDS;
import hdf.object.h5.H5Datatype;
import hdf.object.h5.H5Group;
import hdf.object.h5.H5ScalarDS;
import hdf.object.h5.H5File;
import svljava.SVLJava;
import svljava.SVLJavaException;
import svljava.SVLVar;
import svljava.SVLWriter;
import com.quantumbioinc.datacorrespondent.Correspondent;
//import com.quantumbioinc.operation.TopologySymmetry;
import com.quantumbioinc.xml.Element;
import com.quantumbioinc.xml.Hamiltonian;
import com.quantumbioinc.xml.bind.marshaller.DivconNamespacePrefixMapper;
import com.quantumbioinc.xml.divcon.Atom;
import com.quantumbioinc.xml.divcon.AtomArray;
import com.quantumbioinc.xml.divcon.AtomType;
import com.quantumbioinc.xml.divcon.ChargesType;
import com.quantumbioinc.xml.divcon.Cml;
import com.quantumbioinc.xml.divcon.DivconType;
import com.quantumbioinc.xml.divcon.HamiltonianType;
import com.quantumbioinc.xml.divcon.Molecule;
import com.quantumbioinc.xml.divcon.ObjectFactory;
import com.quantumbioinc.xml.divcon.Scalar;
import com.quantumbioinc.xml.divcon.TotalChargeType;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.math.BigInteger;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.transform.stream.StreamSource;
import hdf.hdf5lib.exceptions.HDF5Exception;
import hdf.object.Group;

        
/**
 *
 * @author roger
 */
public class HDF5Correspondent extends Correspondent implements SVLJavaDriver {

         static {
//         System.out.println(System.getProperty("java.library.path"));
//@todo         System.loadLibrary("JNIDivcon");
//         System.load("/home/roger/NetBeansProjects/JNIDivcon/native/JNIDivcon/dist/Release/GNU-Linux-x86/libJNIDivcon.so");
         }
    private int taskID = 0;		// SVL task identifier
    java.io.PrintStream sessionOutBuffer=null;
    java.io.PrintStream sessionErrBuffer=null;
    
    @Override
    public SVLVar run(SVLVar arg, int tid) throws SVLJavaException
    {
	taskID = tid;

	try {
	    String cmd = arg.getTokn(1);	// subcommand
            SVLJava.print("cmd="+cmd);
	    SVLVar res = arg.peek(2);		// its argument
            SVLJava.print("length="+res.length());
            switch(cmd)
            {
                case "listModels":
                res = listModels(res);
                    break;
                case "storeModel":
                 res = storeModel(res);
                    break;
                case "retrieveModel":
                res = retrieveModel(res);
                    break;
                case "retrievePosingModel":
                res = retrievePosingModel(res);
                    break;
                case "retrieveScalars":
                    res=retrieveScalars(res);
                    break;
                case "retrieveQMScore":
                res = retrieveQMScore(res);
                    break;
                case "retrieveResiduePWD":
                res = retrieveResiduePWD(res);
                    break;
                case "retrieveAtomByAtomPWD":
                res = retrieveAtomByAtomPWD(res);
                    break;
                case "retrieveAtomByAtomDecomposition":
                res = retrieveAtomByAtomDecomposition(res);
                    break;
                case "retrieveAtomByAtomMRM":
                res = retrieveAtomByAtomMRM(res);
                    break;
                case "retrieveNMRScore":
                res = retrieveNMRScore(res);
                    break;
                case "retrieveChemicalShifts":
                res = retrieveChemicalShifts(res);
                    break;
                case "retrieveDensities":
                res = retrieveDensities(res);
                    break;
                case "retrieveEigenVectors":
                res = retrieveEigenVectors(res);
                    break;
                case "retrieveEnergyLevels":
                res = retrieveEnergyLevels(res);
                    break;
                case "setNMRAtomSelection":
                res = setNMRAtomSelection(res);
                    break;
                case "retrieveHamiltonian":
                res = retrieveHamiltonian(res);
                    break;
                case "retrieveDefaultProgramOptions":
                res = retrieveDefaultProgramOptions(res);
                    break;
                case "retrieveLigandSelection":
                res = retrieveLigandSelection(res);
                    break;
                case "retrieveTopologyClassNumbers":
                res = retrieveTopologyClassNumbers(res);
                    break;
                default:
                  res=new SVLVar(new SVLVar(), new SVLVar("Exception:\n"+"here Unknown command: '" + cmd + "'."));
		//throw new SVLJavaException("here Unknown command: '" + cmd + "'.");
                    break;
            }
	    return res;
	}
	catch (Exception ex) {
	    StringWriter sw = new StringWriter();
	    ex.printStackTrace(new PrintWriter(sw));
            return new SVLVar(new SVLVar(), new SVLVar("Exception:\n"+sw.toString()));
//	    throw new SVLJavaException("Exception:\n"+sw.toString());
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
            String nmrXPath="/DivCon/"+targetList.get(index)+"/NMR Score Atom Specific";
            HObject nmrHO=findHDF5Object(h5File, nmrXPath);
            String sdfName="";
            H5Group documentGroup=(H5Group)h5File.get("DivCon/"+targetList.get(index)+"/Documents");
            if(documentGroup!=null)
            {
                ListIterator<HObject> li=documentGroup.getMemberList().listIterator();
                boolean hit=false;
                while(!hit && li.hasNext())
                {
                    HObject hobj=li.next();
                    if(hobj.getName().endsWith(".sdf"))
                    {
                        sdfName=hobj.getName().substring(0, hobj.getName().length()-4);
                        hit=true;
                    }
                }
            }
            H5CompoundDS nmrHObject=(H5CompoundDS)nmrHO;
            if(hObject!=null){
                Vector o=(Vector)hObject.getData();
                individuals[1]=new SVLVar((String[])o.elementAt(0));
            }
            else if(nmrHObject!=null){
                Vector o=(Vector)nmrHObject.getData();
                individuals[1]=new SVLVar((String[])o.elementAt(0));
            }
            else if(!sdfName.isEmpty())
            {
                HObject obj=h5File.get("DivCon/"+targetList.get(index)+"/Documents/"+sdfName+".xml");
        //        java.lang.System.out.println("looking for: "+"DivCon/"+target+"/Documents/"+ligand+".xml");
                if(obj!=null)
                {
                    H5ScalarDS doc=(H5ScalarDS)obj;
                    JAXBContext jc = JAXBContext.newInstance("com.quantumbioinc.xml.divcon");
                    Unmarshaller um=jc.createUnmarshaller();
                    JAXBElement<com.quantumbioinc.xml.divcon.DivconType> jaxbOutElement=um.unmarshal(new StreamSource(new ByteArrayInputStream (((String[])doc.read())[0].getBytes())), com.quantumbioinc.xml.divcon.DivconType.class);
                    DivconType divcon=jaxbOutElement.getValue();
                    Cml cml=(Cml)divcon.getCmlOrTargetOrLigand().get(0);
                    ArrayList<String> titleList=new ArrayList<String>();
                    String[] titles=new String[cml.getAnyCmlOrAnyOrAny().size()];
                    for(int moleculeCount=0;moleculeCount<cml.getAnyCmlOrAnyOrAny().size();moleculeCount++)
                    {
                        JAXBElement<com.quantumbioinc.xml.divcon.Molecule> jaxbMoleculeElement=(JAXBElement<com.quantumbioinc.xml.divcon.Molecule>)cml.getAnyCmlOrAnyOrAny().get(moleculeCount);
                        Molecule molecule=jaxbMoleculeElement.getValue();
                        if(molecule.getTitle()!=null)
                        {
                            if(!titleList.contains(molecule.getTitle()))
                            {
                            titleList.add(molecule.getTitle());
                            }
                            else
                            {
                                titleList.add(molecule.getTitle()+"."+(moleculeCount+1));
                            }
                        }
                        else
                        {
                            titleList.add("Ligand"+(moleculeCount+1));
                        }
                    }
                    for(int moleculeCount=0;moleculeCount<cml.getAnyCmlOrAnyOrAny().size();moleculeCount++)
                    {
                        titles[moleculeCount]=titleList.get(moleculeCount);
                    }
                    individuals[1]=new SVLVar(titles);
                }
                else
                {
                    individuals[1]=new SVLVar();
                }
            }
            else
            {
                individuals[1]=new SVLVar();
            }
            data[index]=new SVLVar(new String[]{"target", "ligand"}, individuals);
        }
        h5File.close();
        return new SVLVar(new SVLVar(data), new SVLVar("",true));
    }

    private SVLVar storeModel(SVLVar var) throws SVLJavaException, IOException, Exception
    {
        String filename = var.peek(1).getTokn(1);
        String target = var.peek(1).getTokn(2);
                sessionErrBuffer=new PrintStream(new java.io.ByteArrayOutputStream(), true);
                sessionOutBuffer=new PrintStream(new java.io.ByteArrayOutputStream(), true);
        java.lang.System.setErr(sessionErrBuffer);
        java.lang.System.setOut(sessionOutBuffer);
        H5File h5File = new H5File(filename, H5File.FILE_CREATE_OPEN);
        if(h5File.exists())
        {
            h5File=(H5File)h5File.createInstance(filename, H5File.WRITE);
        }
        else
        {
            h5File=(H5File)h5File.createInstance(filename, H5File.CREATE);
        }
        System.out.println("h5File.canWrite="+h5File.canWrite());
        h5File.open();
        long[] dims = {1};
        Group targetGroup=(Group)h5File.get("DivCon/"+target+"/Documents");
        if(targetGroup==null)
        {
            Group divconGroup=(Group)h5File.get("DivCon");
            if(divconGroup==null)
            {
                divconGroup=h5File.createGroup("DivCon", null);
            }
            targetGroup=h5File.createGroup(target, divconGroup);
        }
                    ObjectFactory objectFactory=new ObjectFactory();
                    DivconType divcon=new DivconType();
                    divcon.setVersion("1.0");
                    SVLVar svlMol=var.peek(1).peek(3);
                    String[] chainNames=svlMol.peek(2).peek(1).getTokns();
                    int[] residueCounts=svlMol.peek(2).peek(4).getInts();
                    String[] residueNames=svlMol.peek(3).peek(1).getTokns();
                    int[] sequences=svlMol.peek(3).peek(2).getInts();
                    int[] atomCounts=svlMol.peek(3).peek(5).getInts();
                    String[] symbols=svlMol.peek(4).peek(1).getTokns();
                    String[] hybridizations=svlMol.peek(4).peek(3).getTokns();
//                    int[] atomCounts=svlMol.peek(4).peek(5).getInts();
                    String[] atomNames=svlMol.peek(4).peek(8).getTokns();
                    int[] formalCharges=svlMol.peek(4).peek(2).getInts();
                    double[] x=svlMol.peek(4).peek(10).getReals();
                    double[] y=svlMol.peek(4).peek(11).getReals();
                    double[] z=svlMol.peek(4).peek(12).getReals();
                    Cml cml=objectFactory.createCml();
                    int chainIndex=0;
                    int residueIndex=0;
                    int chainOffset=0;
                    int residueOffset=0;
                    JAXBElement<Molecule> molecule=objectFactory.createMolecule(new Molecule());
                    JAXBElement<Molecule> submolecule=objectFactory.createMolecule(new Molecule());
                    JAXBElement<AtomArray> atomArray=objectFactory.createAtomArray(new AtomArray());
                    int totalCharge=0;
                    for(int atomIndex=0;atomIndex<atomNames.length;atomIndex++)
                    {
                        JAXBElement<Atom> atom=objectFactory.createAtom(new Atom());
                        atom.getValue().setElementType(symbols[atomIndex]);
                        atom.getValue().setX3(new Double(x[atomIndex]));
                        atom.getValue().setY3(new Double(y[atomIndex]));
                        atom.getValue().setZ3(new Double(z[atomIndex]));
                        atom.getValue().setFormalCharge(new BigInteger(""+formalCharges[atomIndex]));
                        totalCharge+=formalCharges[atomIndex];
                        JAXBElement<AtomType> atomType=objectFactory.createAtomType(new AtomType());
                        atomType.getValue().setName(atomNames[atomIndex]);
                        atom.getValue().getAnyCmlOrAnyOrAny().add(atomType);
                        JAXBElement<Scalar> chainID=objectFactory.createScalar(new Scalar());
                        chainID.getValue().setTitle("chainID");
                        chainID.getValue().setValue(chainNames[chainIndex].substring(chainNames[chainIndex].lastIndexOf('.')+1, chainNames[chainIndex].length()));
                        atom.getValue().getAnyCmlOrAnyOrAny().add(chainID);
                        JAXBElement<Scalar> hybridization=objectFactory.createScalar(new Scalar());
                        hybridization.getValue().setTitle("hybridization");
                        hybridization.getValue().setValue(hybridizations[atomIndex]);
                        atom.getValue().getAnyCmlOrAnyOrAny().add(hybridization);
                        atomArray.getValue().getAnyCmlOrAnyOrAny().add(atom);
                        if(atomIndex>=atomCounts[residueIndex]+residueOffset-1)
                        {
                            residueOffset+=atomCounts[residueIndex];
                            submolecule.getValue().setTitle(residueNames[residueIndex]);
                            submolecule.getValue().getAnyCmlOrAnyOrAny().add(atomArray);
                            JAXBElement<Scalar> sequence=objectFactory.createScalar(new Scalar());
                            sequence.getValue().setTitle("sequence");
                            sequence.getValue().setValue(""+sequences[residueIndex]);
                            submolecule.getValue().getAnyCmlOrAnyOrAny().add(sequence);
                            molecule.getValue().getAnyCmlOrAnyOrAny().add(submolecule);
                            residueIndex++;
                            if(residueIndex<=atomCounts.length)
                            {
                                atomArray=objectFactory.createAtomArray(new AtomArray());
                                submolecule=objectFactory.createMolecule(new Molecule());
                            }
                            if(residueIndex>=residueCounts[chainIndex]+chainOffset)
                            {
                                chainOffset+=residueCounts[chainIndex];
                                chainIndex++;
                                molecule.getValue().setTitle(svlMol.peek(1).getTokn(1));
                                cml.getAnyCmlOrAnyOrAny().add(molecule);
                                if(chainIndex<=residueCounts.length)molecule=objectFactory.createMolecule(new Molecule());
                            }
                        }
                    }
                    HamiltonianType hamiltonianType=objectFactory.createHamiltonianType();
                    hamiltonianType.setParameters("pm6");
                    divcon.setHamiltonian(hamiltonianType);
                    TotalChargeType totalChargeType=objectFactory.createTotalChargeType();
                    totalChargeType.setValue(totalCharge);
                    ChargesType chargesType=objectFactory.createChargesType();
                    chargesType.setTotalCharge(totalChargeType);
                    divcon.setCharges(chargesType);
                    divcon.getCmlOrTargetOrLigand().add(cml);
                    JAXBContext jc = JAXBContext.newInstance("com.quantumbioinc.xml.divcon");
                    Marshaller m = jc.createMarshaller();
                    m.setProperty(Marshaller.JAXB_SCHEMA_LOCATION, "http://quantumbioinc.com/schema/divcon /home/roger/NetBeansProjects/OOBackbone/schemas/divcon.xsd");
                    m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
                    DivconNamespacePrefixMapper dnpm=new DivconNamespacePrefixMapper();
                    m.setProperty("com.sun.xml.bind.namespacePrefixMapper", dnpm);
            //Marshal object into file.
                    ByteArrayOutputStream baos=new java.io.ByteArrayOutputStream();
                    JAXBElement<com.quantumbioinc.xml.divcon.DivconType> jaxbOutElement=objectFactory.createDivcon(divcon);
                    m.marshal(jaxbOutElement, baos);
        Datatype docType=h5File.createDatatype(H5Datatype.CLASS_STRING, baos.size(), Datatype.NATIVE, Datatype.NATIVE);
        ListIterator<HObject> li=targetGroup.getMemberList().listIterator();
        boolean documentExists=false;
        while(li.hasNext())
        {
            Object obj=li.next();
            if(obj instanceof hdf.object.h5.H5ScalarDS && ((hdf.object.h5.H5ScalarDS)obj).getName().compareTo(target+".xml")==0)
            {
                documentExists=true;
                ((hdf.object.h5.H5ScalarDS)obj).write(new String[]{baos.toString()});
                break;
            }
        }
        if(!documentExists)h5File.createScalarDS(target+".xml", targetGroup, docType, dims, null, null, 0, new String[]{baos.toString()});
        h5File.close();
        return new SVLVar();
    }

    private SVLVar retrieveModel(SVLVar var) throws SVLJavaException, IOException, Exception
    {
        String filename = var.peek(1).getTokn(1);
        String target = var.peek(1).getTokn(2);
                sessionErrBuffer=new PrintStream(new java.io.ByteArrayOutputStream(), true);
                sessionOutBuffer=new PrintStream(new java.io.ByteArrayOutputStream(), true);
        java.lang.System.setErr(sessionErrBuffer);
        java.lang.System.setOut(sessionOutBuffer);
//        java.lang.System.out.println("in: retrieveModel");
        H5File h5File = new H5File(filename, H5File.READ);
        if(!h5File.exists())
        {
            return new SVLVar(new SVLVar(), new SVLVar(filename + " doesn't exists.", true));
//            throw new SVLJavaException(filename + " doesn't exists.");
        }
        h5File=(H5File)h5File.createInstance(filename, H5File.READ);
        h5File.open();
        long[] dims = {1};
        Group targetGroup=(Group)h5File.get("DivCon/"+target);
        if(targetGroup==null)
        {
            return new SVLVar(new SVLVar(), new SVLVar(target + " specimen doesn't exists in "+filename+".", true));
            //throw new SVLJavaException(target + " specimen doesn't exists in "+filename+".");
        }
        HObject obj=h5File.get("DivCon/"+target+"/Documents/"+target+".xml");
//        java.lang.System.out.println("looking for: "+"DivCon/"+target+"/Documents/"+target+".xml");
        SVLVar[] model=new SVLVar[4];
        if(obj!=null)
        {
            H5ScalarDS doc=(H5ScalarDS)obj;
        JAXBContext jc = JAXBContext.newInstance("com.quantumbioinc.xml.divcon");
        Unmarshaller um=jc.createUnmarshaller();
        JAXBElement<com.quantumbioinc.xml.divcon.DivconType> jaxbOutElement=um.unmarshal(new StreamSource(new ByteArrayInputStream (((String[])doc.read())[0].getBytes())), com.quantumbioinc.xml.divcon.DivconType.class);
        DivconType divcon=jaxbOutElement.getValue();
        Cml cml=(Cml)divcon.getCmlOrTargetOrLigand().get(0);
        ArrayList<String> chainTitlesList=new ArrayList<>();
        ArrayList<String> chainsList=new ArrayList<>();
        ArrayList<Integer> chainResidueCountsList=new ArrayList<>();
        ArrayList<String> residueNamesList=new ArrayList<>();
        ArrayList<String> residueAminosList=new ArrayList<>();
        ArrayList<String> aminosList=new ArrayList<>();
        aminosList.add("ALA");
        aminosList.add("ARG");
        aminosList.add("ASN");
        aminosList.add("ALP");
        aminosList.add("CYS");
        aminosList.add("GLN");
        aminosList.add("GLU");
        aminosList.add("GLY");
        aminosList.add("HIS");
        aminosList.add("ILE");
        aminosList.add("LEU");
        aminosList.add("LYS");
        aminosList.add("MET");
        aminosList.add("PHE");
        aminosList.add("PRO");
        aminosList.add("SER");
        aminosList.add("THR");
        aminosList.add("TRP");
        aminosList.add("VAL");
        ArrayList<Integer> residueAtomCountsList=new ArrayList<>();
        ArrayList<Integer> residueOffsetsList=new ArrayList<>();
        ArrayList<Integer> sequencesList=new ArrayList<>();
        String insertionCode="";
        ArrayList<String> namesList=new ArrayList<>();
        ArrayList<String> symbolsList=new ArrayList<>();
        ArrayList<Double> xsList=new ArrayList<>();
        ArrayList<Double> ysList=new ArrayList<>();
        ArrayList<Double> zsList=new ArrayList<>();
        ArrayList<Double> formalChargesList=new ArrayList<>();
        ArrayList<String> hybridizationsList=new ArrayList<>();
        String title=target;
        for(int moleculeCount=0;moleculeCount<cml.getAnyCmlOrAnyOrAny().size();moleculeCount++)
        {
            if(((JAXBElement)cml.getAnyCmlOrAnyOrAny().get(moleculeCount)).getValue() instanceof com.quantumbioinc.xml.divcon.Symmetry) continue;
            JAXBElement<com.quantumbioinc.xml.divcon.Molecule> jaxbMoleculeElement=(JAXBElement<com.quantumbioinc.xml.divcon.Molecule>)cml.getAnyCmlOrAnyOrAny().get(moleculeCount);
            Molecule molecule=jaxbMoleculeElement.getValue();
            if(molecule.getTitle()!=null)title=molecule.getTitle();
            String chain="";
            int residueCounts=molecule.getAnyCmlOrAnyOrAny().size();
            for(int submoleculeCount=0;submoleculeCount<molecule.getAnyCmlOrAnyOrAny().size();submoleculeCount++)
            {
                if(((JAXBElement)molecule.getAnyCmlOrAnyOrAny().get(submoleculeCount)).getValue() instanceof com.quantumbioinc.xml.divcon.BondArray){residueCounts--;continue;};
                JAXBElement<com.quantumbioinc.xml.divcon.Molecule> jaxbSubmoleculeElement=(JAXBElement<com.quantumbioinc.xml.divcon.Molecule>)molecule.getAnyCmlOrAnyOrAny().get(submoleculeCount);
                Molecule submolecule=jaxbSubmoleculeElement.getValue();
                JAXBElement<com.quantumbioinc.xml.divcon.AtomArray> jaxbAtomArrayElement=(JAXBElement<com.quantumbioinc.xml.divcon.AtomArray>)submolecule.getAnyCmlOrAnyOrAny().get(0);
                AtomArray atomArray=jaxbAtomArrayElement.getValue();
                for(int atomCount=0;atomCount<atomArray.getAnyCmlOrAnyOrAny().size();atomCount++)
                {
                    JAXBElement<com.quantumbioinc.xml.divcon.Atom> jaxbAtomElement=(JAXBElement<com.quantumbioinc.xml.divcon.Atom>)atomArray.getAnyCmlOrAnyOrAny().get(atomCount);
                    Atom atom=jaxbAtomElement.getValue();
                    if(atom.getAnyCmlOrAnyOrAny().size()>0 && ((JAXBElement)atom.getAnyCmlOrAnyOrAny().get(0)).getValue() instanceof com.quantumbioinc.xml.divcon.AtomType)
                    {
                        JAXBElement<com.quantumbioinc.xml.divcon.AtomType> jaxbAtomTypeElement=(JAXBElement<com.quantumbioinc.xml.divcon.AtomType>)atom.getAnyCmlOrAnyOrAny().get(0);
                        namesList.add(jaxbAtomTypeElement.getValue().getName());
                    }
                    else
                    {
                        namesList.add(atom.getElementType());
                    }
                    symbolsList.add(atom.getElementType());
                    xsList.add(new Double(atom.getX3().doubleValue()));
                    ysList.add(new Double(atom.getY3().doubleValue()));
                    zsList.add(new Double(atom.getZ3().doubleValue()));
                    if(atom.getFormalCharge()!=null)
                    {
                        formalChargesList.add(new Double(atom.getFormalCharge().doubleValue()));
                    }
                    else
                    {
                        formalChargesList.add(new Double(0.0));
                    }
                    if(atom.getAnyCmlOrAnyOrAny().size()>0)
                    {
                        boolean hit=false;
                        int count=0;
                        while(!hit && count<atom.getAnyCmlOrAnyOrAny().size())
                        {
                            if(((JAXBElement)atom.getAnyCmlOrAnyOrAny().get(count)).getValue() instanceof com.quantumbioinc.xml.divcon.Scalar)
                            {
                        JAXBElement<com.quantumbioinc.xml.divcon.Scalar> jaxScalarElement=(JAXBElement<com.quantumbioinc.xml.divcon.Scalar>)atom.getAnyCmlOrAnyOrAny().get(count);
                        switch (jaxScalarElement.getValue().getTitle())
                        {
                            case "chainID":
                                chain=(jaxScalarElement.getValue().getValue());
                                hit=true;
                                break;
                        }
                            }
                            count++;
                        }
                        if(!hit)chain+="";
                    }
                    if(atom.getAnyCmlOrAnyOrAny().size()>0)
                    {
                        boolean hit=false;
                        int count=0;
                        while(!hit && count<atom.getAnyCmlOrAnyOrAny().size())
                        {
                            if(((JAXBElement)atom.getAnyCmlOrAnyOrAny().get(count)).getValue() instanceof com.quantumbioinc.xml.divcon.Scalar)
                            {
                        JAXBElement<com.quantumbioinc.xml.divcon.Scalar> jaxScalarElement=(JAXBElement<com.quantumbioinc.xml.divcon.Scalar>)atom.getAnyCmlOrAnyOrAny().get(count);
                        switch (jaxScalarElement.getValue().getTitle())
                        {
                            case "hybridization":
                        hybridizationsList.add(jaxScalarElement.getValue().getValue());
                        hit=true;
                                break;
                        }
                            }
                            count++;
                        }
                        if(!hit)hybridizationsList.add("huh");
                    }
                    if(atomCount==0 && atom.getAnyCmlOrAnyOrAny().size()>0)
                    {
                        boolean hit=false;
                        int count=0;
                        while(!hit && count<atom.getAnyCmlOrAnyOrAny().size())
                        {
                            if(((JAXBElement)atom.getAnyCmlOrAnyOrAny().get(count)).getValue() instanceof com.quantumbioinc.xml.divcon.Scalar)
                            {
                        JAXBElement<com.quantumbioinc.xml.divcon.Scalar> jaxScalarElement=(JAXBElement<com.quantumbioinc.xml.divcon.Scalar>)atom.getAnyCmlOrAnyOrAny().get(count);
                        switch (jaxScalarElement.getValue().getTitle())
                        {
                            case "insertionCode":
                                insertionCode+=jaxScalarElement.getValue().getValue();
                                hit=true;
                                break;
                        }
                            }
                        count++;
                        }
                        if(!hit)insertionCode+=" ";
                    }
                    else if(atomCount==0)
                    {
                        insertionCode+=" ";
                    }
                }
                if(submolecule.getAnyCmlOrAnyOrAny().size()>0)
                {
                        boolean hit=false;
                        int count=1;
                        while(!hit && count<submolecule.getAnyCmlOrAnyOrAny().size())
                        {
                            if(((JAXBElement)submolecule.getAnyCmlOrAnyOrAny().get(count)).getValue() instanceof com.quantumbioinc.xml.divcon.Scalar)
                            {
                                JAXBElement<com.quantumbioinc.xml.divcon.Scalar> jaxScalarElement=(JAXBElement<com.quantumbioinc.xml.divcon.Scalar>)submolecule.getAnyCmlOrAnyOrAny().get(count);
                                switch (jaxScalarElement.getValue().getTitle())
                                {
                                    case "sequence":
                                        sequencesList.add(new Integer(jaxScalarElement.getValue().getValue()));
                                        hit=true;
                                        break;
                                }
                            }
                            count++;
                        }
                        if(!hit)sequencesList.add(new Integer(0));
                }
                if(submolecule.getTitle()!=null)
                {
                    residueNamesList.add(submolecule.getTitle());
                }
                else
                {
                    residueNamesList.add("");
                }
                residueAtomCountsList.add(new Integer(atomArray.getAnyCmlOrAnyOrAny().size()));
                if(aminosList.contains(submolecule.getTitle()))
                {
                    residueAminosList.add("amino");
                }
                else
                {
                    residueAminosList.add("none");
                }
            }
            chainTitlesList.add(title);
            chainsList.add(title+"."+chain);
            chainResidueCountsList.add(new Integer(residueCounts));
        }
        String[] jchainTitles=new String[chainTitlesList.size()];
        chainTitlesList.toArray(jchainTitles);
        String[] jchains=new String[chainsList.size()];
        chainsList.toArray(jchains);
        String[] jchainEmpties=new String[chainsList.size()];
        int[] jchainResidueCounts=new int[chainsList.size()];
        for(int vIndex=0;vIndex<jchainResidueCounts.length;vIndex++)
        {
            jchainEmpties[vIndex]="";
            jchainResidueCounts[vIndex]=chainResidueCountsList.get(vIndex);
        }
        String[] jresidueNames=new String[residueNamesList.size()];
        residueNamesList.toArray(jresidueNames);
        String[] jresidueAminos=new String[residueAminosList.size()];
        residueAminosList.toArray(jresidueAminos);
        int[] jresidueAtomCounts=new int[residueAtomCountsList.size()];
        int[] jsequences=new int[residueAtomCountsList.size()];
        for(int vIndex=0;vIndex<jresidueAtomCounts.length;vIndex++)
        {
            jsequences[vIndex]=sequencesList.get(vIndex);
            jresidueAtomCounts[vIndex]=residueAtomCountsList.get(vIndex);
        }
        String[] jnames=new String[namesList.size()];
        namesList.toArray(jnames);
        String[] jsymbols=new String[symbolsList.size()];
        symbolsList.toArray(jsymbols);
        String[] jhybridizations=new String[hybridizationsList.size()];
        hybridizationsList.toArray(jhybridizations);
        SVLVar[] linkageVectors=new SVLVar[4];
        linkageVectors[0]=new SVLVar(jchains);
        linkageVectors[1]=new SVLVar(jchainTitles);
        linkageVectors[2]=new SVLVar(jchainEmpties);
        linkageVectors[3]=new SVLVar(jchainResidueCounts);
        SVLVar[] residueVectors=new SVLVar[5];
        residueVectors[0]=new SVLVar(jresidueNames);
        residueVectors[1]=new SVLVar(jsequences);
        residueVectors[2]=new SVLVar(insertionCode);
        residueVectors[3]=new SVLVar(jresidueAminos);
        residueVectors[4]=new SVLVar(jresidueAtomCounts);
        //if(true) return new SVLVar(new SVLVar(residueVectors[0]), new SVLVar("",true));
        SVLVar[] atomVectors=new SVLVar[12];
        double[] jxs=new double[xsList.size()];
        double[] jys=new double[ysList.size()];
        double[] jzs=new double[zsList.size()];
        double[] jformalCharges=new double[formalChargesList.size()];
        for(int vIndex=0;vIndex<jformalCharges.length;vIndex++)
        {
            jxs[vIndex]=xsList.get(vIndex);
            jys[vIndex]=ysList.get(vIndex);
            jzs[vIndex]=zsList.get(vIndex);
            jformalCharges[vIndex]=formalChargesList.get(vIndex);
        }
        atomVectors[0]=new SVLVar(jsymbols);
        atomVectors[1]=new SVLVar(jformalCharges);
        atomVectors[2]=new SVLVar(jhybridizations);
        atomVectors[3]=new SVLVar("");
        atomVectors[4]=new SVLVar("");
        atomVectors[5]=new SVLVar("");
        atomVectors[6]=new SVLVar("");
        atomVectors[7]=new SVLVar(jnames);
        atomVectors[8]=new SVLVar("");
        atomVectors[9]=new SVLVar(jxs);
        atomVectors[10]=new SVLVar(jys);
        atomVectors[11]=new SVLVar(jzs);
        model[0]=new SVLVar(title, true);
        model[1]=new SVLVar(linkageVectors);
        model[2]=new SVLVar(residueVectors);
        model[3]=new SVLVar(atomVectors);
        }
        h5File.close();
        return new SVLVar(new SVLVar(new SVLVar(model)), new SVLVar("",true));
    }

    private SVLVar retrievePosingModel(SVLVar var) throws SVLJavaException, IOException, Exception
    {
        String filename = var.peek(1).getTokn(1);
        String target = var.peek(1).getTokn(2);
        String ligand = var.peek(1).getTokn(3);
                sessionErrBuffer=new PrintStream(new java.io.ByteArrayOutputStream(), true);
                sessionOutBuffer=new PrintStream(new java.io.ByteArrayOutputStream(), true);
        java.lang.System.setErr(sessionErrBuffer);
        java.lang.System.setOut(sessionOutBuffer);
//        java.lang.System.out.println("in: retrieveModel");
        H5File h5File = new H5File(filename, H5File.READ);
        if(!h5File.exists())
        {
            return new SVLVar(new SVLVar(), new SVLVar(filename + " doesn't exists.", true));
//            throw new SVLJavaException(filename + " doesn't exists.");
        }
        h5File=(H5File)h5File.createInstance(filename, H5File.READ);
        h5File.open();
        long[] dims = {1};
        Group targetGroup=(Group)h5File.get("DivCon/"+target);
        if(targetGroup==null)
        {
            return new SVLVar(new SVLVar(), new SVLVar(target + " specimen doesn't exists in "+filename+".", true));
            //throw new SVLJavaException(target + " specimen doesn't exists in "+filename+".");
        }
        String sdfName="";
        H5Group documentGroup=(H5Group)h5File.get("DivCon/"+target+"/Documents");
        if(documentGroup!=null)
        {
            ListIterator<HObject> li=documentGroup.getMemberList().listIterator();
            boolean hit=false;
            while(!hit && li.hasNext())
            {
                HObject ho=li.next();
                if(ho.getName().endsWith(".sdf"))
                {
                    sdfName=ho.getName().substring(0, ho.getName().length()-4);
                    hit=true;
                }
            }
            if(!hit)
            {
                sdfName=ligand;
            }
        }
        HObject obj=h5File.get("DivCon/"+target+"/Documents/"+sdfName+".xml");
//        java.lang.System.out.println("looking for: "+"DivCon/"+target+"/Documents/"+ligand+".xml");
        SVLVar[] model=new SVLVar[4];
        if(obj!=null)
        {
            H5ScalarDS doc=(H5ScalarDS)obj;
        JAXBContext jc = JAXBContext.newInstance("com.quantumbioinc.xml.divcon");
        Unmarshaller um=jc.createUnmarshaller();
        JAXBElement<com.quantumbioinc.xml.divcon.DivconType> jaxbOutElement=um.unmarshal(new StreamSource(new ByteArrayInputStream (((String[])doc.read())[0].getBytes())), com.quantumbioinc.xml.divcon.DivconType.class);
        DivconType divcon=jaxbOutElement.getValue();
        Cml cml=(Cml)divcon.getCmlOrTargetOrLigand().get(0);
        ArrayList<String> chainTitlesList=new ArrayList<>();
        ArrayList<String> chainsList=new ArrayList<>();
        ArrayList<Integer> chainResidueCountsList=new ArrayList<>();
        ArrayList<String> residueNamesList=new ArrayList<>();
        ArrayList<String> residueAminosList=new ArrayList<>();
        ArrayList<String> aminosList=new ArrayList<>();
        aminosList.add("ALA");
        aminosList.add("ARG");
        aminosList.add("ASN");
        aminosList.add("ALP");
        aminosList.add("CYS");
        aminosList.add("GLN");
        aminosList.add("GLU");
        aminosList.add("GLY");
        aminosList.add("HIS");
        aminosList.add("ILE");
        aminosList.add("LEU");
        aminosList.add("LYS");
        aminosList.add("MET");
        aminosList.add("PHE");
        aminosList.add("PRO");
        aminosList.add("SER");
        aminosList.add("THR");
        aminosList.add("TRP");
        aminosList.add("VAL");
        ArrayList<Integer> residueAtomCountsList=new ArrayList<>();
        ArrayList<Integer> residueOffsetsList=new ArrayList<>();
        ArrayList<Integer> sequencesList=new ArrayList<>();
        String insertionCode="";        
        ArrayList<String> namesList=new ArrayList<>();
        ArrayList<String> symbolsList=new ArrayList<>();
        ArrayList<Double> xsList=new ArrayList<>();
        ArrayList<Double> ysList=new ArrayList<>();
        ArrayList<Double> zsList=new ArrayList<>();
        ArrayList<Double> formalChargesList=new ArrayList<>();
        ArrayList<String> hybridizationsList=new ArrayList<>();
        String title="";
        int moleculeCount=0;
        boolean hitPose=false;
        while(!hitPose && moleculeCount<cml.getAnyCmlOrAnyOrAny().size())
        {
            JAXBElement<com.quantumbioinc.xml.divcon.Molecule> jaxbMoleculeElement=(JAXBElement<com.quantumbioinc.xml.divcon.Molecule>)cml.getAnyCmlOrAnyOrAny().get(moleculeCount);
            Molecule molecule=jaxbMoleculeElement.getValue();
//            if(molecule.getTitle()==null || molecule.getTitle().compareTo(ligand)!=0){moleculeCount++;continue;};
            if(molecule.getTitle()!=null && molecule.getTitle().compareTo(ligand)==0)
            {
                hitPose=true;
                title=molecule.getTitle();
            }
            else
            {
                moleculeCount++;
            }
        }
        if(!hitPose)
        {
            if(ligand.lastIndexOf('.')>=0)
            {
                try  
                {  
                    int d = Integer.parseInt(ligand.substring(ligand.lastIndexOf('.')+1));
                    hitPose=true;
                    moleculeCount=d-1;
                    title=ligand.substring(0, ligand.lastIndexOf('.'));
                }  
                catch(NumberFormatException nfe)  
                {  
                }
            }
            
        }
        if(hitPose)
        {
            JAXBElement<com.quantumbioinc.xml.divcon.Molecule> jaxbMoleculeElement=(JAXBElement<com.quantumbioinc.xml.divcon.Molecule>)cml.getAnyCmlOrAnyOrAny().get(moleculeCount);
            Molecule molecule=jaxbMoleculeElement.getValue();
            String chain="";
            for(int submoleculeCount=0;submoleculeCount<molecule.getAnyCmlOrAnyOrAny().size();submoleculeCount++)
            {
                JAXBElement<com.quantumbioinc.xml.divcon.Molecule> jaxbSubmoleculeElement=(JAXBElement<com.quantumbioinc.xml.divcon.Molecule>)molecule.getAnyCmlOrAnyOrAny().get(submoleculeCount);
                Molecule submolecule=jaxbSubmoleculeElement.getValue();
                JAXBElement<com.quantumbioinc.xml.divcon.AtomArray> jaxbAtomArrayElement=(JAXBElement<com.quantumbioinc.xml.divcon.AtomArray>)submolecule.getAnyCmlOrAnyOrAny().get(0);
                AtomArray atomArray=jaxbAtomArrayElement.getValue();
                for(int atomCount=0;atomCount<atomArray.getAnyCmlOrAnyOrAny().size();atomCount++)
                {
                    JAXBElement<com.quantumbioinc.xml.divcon.Atom> jaxbAtomElement=(JAXBElement<com.quantumbioinc.xml.divcon.Atom>)atomArray.getAnyCmlOrAnyOrAny().get(atomCount);
                    Atom atom=jaxbAtomElement.getValue();
                    if(atom.getAnyCmlOrAnyOrAny().size()>0 && ((JAXBElement)atom.getAnyCmlOrAnyOrAny().get(0)).getValue() instanceof com.quantumbioinc.xml.divcon.AtomType)
                    {
                    JAXBElement<com.quantumbioinc.xml.divcon.AtomType> jaxbAtomTypeElement=(JAXBElement<com.quantumbioinc.xml.divcon.AtomType>)atom.getAnyCmlOrAnyOrAny().get(0);
                    namesList.add(jaxbAtomTypeElement.getValue().getName());
                    }
                    else
                    {
                        namesList.add("");
                    }
                    symbolsList.add(atom.getElementType());
                    xsList.add(new Double(atom.getX3().doubleValue()));
                    ysList.add(new Double(atom.getY3().doubleValue()));
                    zsList.add(new Double(atom.getZ3().doubleValue()));
                    if(atom.getFormalCharge()!=null)
                    {
                        formalChargesList.add(new Double(atom.getFormalCharge().doubleValue()));
                    }
                    else
                    {
                        formalChargesList.add(new Double(0.0));
                    }
                    if(atom.getAnyCmlOrAnyOrAny().size()>1)
                    {
                        JAXBElement<com.quantumbioinc.xml.divcon.Scalar> jaxScalarElement=(JAXBElement<com.quantumbioinc.xml.divcon.Scalar>)atom.getAnyCmlOrAnyOrAny().get(1);
                        switch (jaxScalarElement.getValue().getTitle())
                        {
                            case "chainID":
                                chain=(jaxScalarElement.getValue().getValue());
                                break;
                        }
                    }
                    if(atom.getAnyCmlOrAnyOrAny().size()>2)
                    {
                        boolean hit=false;
                        int count=0;
                        while(!hit && count<atom.getAnyCmlOrAnyOrAny().size())
                        {
                            if(((JAXBElement)atom.getAnyCmlOrAnyOrAny().get(count)).getValue() instanceof com.quantumbioinc.xml.divcon.Scalar)
                            {
                        JAXBElement<com.quantumbioinc.xml.divcon.Scalar> jaxScalarElement=(JAXBElement<com.quantumbioinc.xml.divcon.Scalar>)atom.getAnyCmlOrAnyOrAny().get(count);
                        switch (jaxScalarElement.getValue().getTitle())
                        {
                            case "hybridization":
                        hybridizationsList.add(jaxScalarElement.getValue().getValue());
                        hit=true;
                                break;
                        }
                            }
                            count++;
                        }
                        if(!hit)hybridizationsList.add("huh");
                    }
                    if(atomCount==0 && atom.getAnyCmlOrAnyOrAny().size()>0)
                    {
                        boolean hit=false;
                        int count=0;
                        while(!hit && count<atom.getAnyCmlOrAnyOrAny().size())
                        {
                            if(((JAXBElement)atom.getAnyCmlOrAnyOrAny().get(count)).getValue() instanceof com.quantumbioinc.xml.divcon.Scalar)
                            {
                        JAXBElement<com.quantumbioinc.xml.divcon.Scalar> jaxScalarElement=(JAXBElement<com.quantumbioinc.xml.divcon.Scalar>)atom.getAnyCmlOrAnyOrAny().get(count);
                        switch (jaxScalarElement.getValue().getTitle())
                        {
                            case "insertionCode":
                                insertionCode+=jaxScalarElement.getValue().getValue();
                                hit=true;
                                break;
                        }
                            }
                        count++;
                        }
                        if(!hit)insertionCode+=" ";
                    }
                    else if(atomCount==0)
                    {
                        insertionCode+=" ";
                    }
                }
                if(submolecule.getAnyCmlOrAnyOrAny().size()>0)
                {
                        boolean hit=false;
                        int count=0;
                        while(!hit && count<submolecule.getAnyCmlOrAnyOrAny().size())
                        {
                            if(((JAXBElement)submolecule.getAnyCmlOrAnyOrAny().get(count)).getValue() instanceof com.quantumbioinc.xml.divcon.Scalar)
                            {
                    JAXBElement<com.quantumbioinc.xml.divcon.Scalar> jaxScalarElement=(JAXBElement<com.quantumbioinc.xml.divcon.Scalar>)submolecule.getAnyCmlOrAnyOrAny().get(count);
                    switch (jaxScalarElement.getValue().getTitle())
                    {
                        case "sequence":
                            sequencesList.add(new Integer(jaxScalarElement.getValue().getValue()));
                            hit=true;
                            break;
                    }
                            }
                        count++;
                        }
                        if(!hit)sequencesList.add(new Integer(0));
                }
//        if(true) return new SVLVar(new SVLVar("here", true), new SVLVar("",true));
                if(submolecule.getTitle()!=null)
                {
                    residueNamesList.add(submolecule.getTitle());
                }
                else
                {
                    residueNamesList.add("");
                }
                residueAtomCountsList.add(new Integer(atomArray.getAnyCmlOrAnyOrAny().size()));
                if(aminosList.contains(residueNamesList.get(residueNamesList.size()-1)))
                {
                    residueAminosList.add("amino");
                }
                else
                {
                    residueAminosList.add("none");
                }
            }
            chainTitlesList.add(title);
            chainsList.add(title+"."+chain);
            chainResidueCountsList.add(new Integer(molecule.getAnyCmlOrAnyOrAny().size()));
        }
        //if(true) new SVLVar(new SVLVar("here"), new SVLVar("",true));
        String[] jchainTitles=new String[chainTitlesList.size()];
        chainTitlesList.toArray(jchainTitles);
        String[] jchains=new String[chainsList.size()];
        chainsList.toArray(jchains);
        String[] jchainEmpties=new String[chainsList.size()];
        int[] jchainResidueCounts=new int[chainsList.size()];
        for(int vIndex=0;vIndex<jchainResidueCounts.length;vIndex++)
        {
            jchainEmpties[vIndex]="";
            jchainResidueCounts[vIndex]=chainResidueCountsList.get(vIndex);
        }
        String[] jresidueNames=new String[residueNamesList.size()];
        residueNamesList.toArray(jresidueNames);
        String[] jresidueAminos=new String[residueAminosList.size()];
        residueAminosList.toArray(jresidueAminos);
        int[] jresidueAtomCounts=new int[residueAtomCountsList.size()];
        int[] jsequences=new int[residueAtomCountsList.size()];
        for(int vIndex=0;vIndex<jresidueAtomCounts.length;vIndex++)
        {
            jsequences[vIndex]=sequencesList.get(vIndex);
            jresidueAtomCounts[vIndex]=residueAtomCountsList.get(vIndex);
        }
        String[] jnames=new String[namesList.size()];
        namesList.toArray(jnames);
        String[] jsymbols=new String[symbolsList.size()];
        symbolsList.toArray(jsymbols);
        String[] jhybridizations=new String[hybridizationsList.size()];
        hybridizationsList.toArray(jhybridizations);
        SVLVar[] linkageVectors=new SVLVar[4];
        linkageVectors[0]=new SVLVar(jchains);
        linkageVectors[1]=new SVLVar(jchainTitles);
        linkageVectors[2]=new SVLVar(jchainEmpties);
        linkageVectors[3]=new SVLVar(jchainResidueCounts);
        SVLVar[] residueVectors=new SVLVar[5];
        residueVectors[0]=new SVLVar(jresidueNames);
        residueVectors[1]=new SVLVar(jsequences);
        residueVectors[2]=new SVLVar(insertionCode);
        residueVectors[3]=new SVLVar(jresidueAminos);
        residueVectors[4]=new SVLVar(jresidueAtomCounts);
        SVLVar[] atomVectors=new SVLVar[12];
        double[] jxs=new double[xsList.size()];
        double[] jys=new double[ysList.size()];
        double[] jzs=new double[zsList.size()];
        double[] jformalCharges=new double[formalChargesList.size()];
        for(int vIndex=0;vIndex<jformalCharges.length;vIndex++)
        {
            jxs[vIndex]=xsList.get(vIndex);
            jys[vIndex]=ysList.get(vIndex);
            jzs[vIndex]=zsList.get(vIndex);
            jformalCharges[vIndex]=formalChargesList.get(vIndex);
        }
        atomVectors[0]=new SVLVar(jsymbols);
        atomVectors[1]=new SVLVar(jformalCharges);
        atomVectors[2]=new SVLVar(jhybridizations);
        atomVectors[3]=new SVLVar("");
        atomVectors[4]=new SVLVar("");
        atomVectors[5]=new SVLVar("");
        atomVectors[6]=new SVLVar("");
        atomVectors[7]=new SVLVar(jnames);
        atomVectors[8]=new SVLVar("");
        atomVectors[9]=new SVLVar(jxs);
        atomVectors[10]=new SVLVar(jys);
        atomVectors[11]=new SVLVar(jzs);
        model[0]=new SVLVar(title, true);
        model[1]=new SVLVar(linkageVectors);
        model[2]=new SVLVar(residueVectors);
        model[3]=new SVLVar(atomVectors);
        }
        h5File.close();
        return new SVLVar(new SVLVar(new SVLVar(model)), new SVLVar("",true));
    }

private SVLVar retrieveScalars(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
//   java.lang.System.out.println("h5 opened ");
    String xPath="/DivCon/"+target;
    H5Group targetGroup=(H5Group)findHDF5Object(h5File, xPath);
    SVLVar[] scalars=new SVLVar[1];
    String[] tags=new String[1];
    if(targetGroup!=null)
    {
        List scalarRow=targetGroup.getMetadata();
        scalars=new SVLVar[scalarRow.size()];
        tags=new String[scalarRow.size()];
        for(int memberIndex=0;memberIndex<scalarRow.size();memberIndex++)
        {
            hdf.object.Attribute member=(hdf.object.Attribute)scalarRow.get(memberIndex);
            tags[memberIndex]=member.getName();
            System.out.println(member.getName()+" value class "+member.getValue().getClass().getName());
            scalars[memberIndex]=new SVLVar(((double[])member.getValue())[0]);
        }
    }
    else
    {
        scalars[0]=new SVLVar(new String[]{filename, target});
    }
    SVLVar data=new SVLVar(tags, scalars);
    h5File.close();
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
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
        if(hObject==null)
        {
            return new SVLVar(new SVLVar(), new SVLVar("QM Score does not exist for: '" + target + "'.", true));
            //throw new SVLJavaException("QM Score does not exist for: '" + target + "'.");
        }
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
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
}

private SVLVar retrieveResiduePWD(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    String ligand="";
    if(var.peek(1).length()>2)
    {
        ligand = var.peek(1).getTokn(3);
        //if(true) return new SVLVar("retrieveAtomByAtomPWD "+ligand);
    }
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
        String xPath="/DivCon/"+target+"/PWD Solvation";
        HObject ho=findHDF5Object(h5File, xPath);
        H5Group hObject=(H5Group)ho;
        if(hObject==null) return new SVLVar();
    if(ligand.length()>0)
    {
            List<HObject> ligandList=hObject.getMemberList();
            H5Group pwdTargetLigandObject=null;
            boolean hit=false;
            for(int index=0;index<hObject.getNumberOfMembersInFile();index++)
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
            SVLVar[] pwdDataset=new SVLVar[pwdTargetLigandObject.getNumberOfMembersInFile()+1];
            if(pwdDataset.length!=5)
            {
                return new SVLVar(new SVLVar(), new SVLVar("pwd set wrong size: " + pwdDataset.length + ".", true));
                //throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            }
            List pwdRow=pwdTargetLigandObject.getMemberList();
            pwdDataset[0]=new SVLVar(pwdTargetLigandObject.getName());
            if(pwdRow.size()==4)
            {
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
            }
            else
            {
                pwdDataset[1]=new SVLVar();
                pwdDataset[2]=new SVLVar();
                pwdDataset[3]=new SVLVar();
                pwdDataset[4]=new SVLVar();
            }
//            pwdDataset[1]=new SVLVar((int[])sequenceA.read());
//            H5ScalarDS sequenceB=(H5ScalarDS)pwdRow.get(1);
//            pwdDataset[2]=new SVLVar((int[])sequenceB.read());
//            H5ScalarDS sequenceValues=(H5ScalarDS)pwdRow.get(2);
//            pwdDataset[3]=new SVLVar((double[])sequenceValues.read());
//            H5ScalarDS sequenceLabels=(H5ScalarDS)pwdRow.get(3);
//            pwdDataset[4]=new SVLVar(sequenceLabels.readBytes());
           SVLVar data=new SVLVar(new String[]{"name", "indexA", "indexB", "values", "labels"}, pwdDataset);
            h5File.close();
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
    }
    else
    {
    SVLVar[] data=new SVLVar[hObject.getNumberOfMembersInFile()];
        List pwdMembers=hObject.getMemberList();
        for(int index=0;index<hObject.getNumberOfMembersInFile();index++)
        {
            H5Group pwdGroup=(H5Group)pwdMembers.get(index);
            SVLVar[] pwdDataset=new SVLVar[pwdGroup.getNumberOfMembersInFile()+1];
            if(pwdDataset.length!=5)
            {
                return new SVLVar(new SVLVar(), new SVLVar("pwd set wrong size: " + pwdDataset.length + "."));
                //throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            }
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
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
    }
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
                        if(obj instanceof hdf.object.Attribute && ((hdf.object.Attribute)obj).getName().compareTo("Total Charge")==0)
                        {
//                            ChargesType chargeType=new ChargesType();
//                            TotalChargeType totalChargeType=new TotalChargeType();
//                            chargeType.setTotalCharge(totalChargeType);
//                            int[] tc=(int[])((hdf.object.Attribute)obj).getValue();
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
            return new SVLVar(new SVLVar(), new SVLVar("", true));
            }
            SVLVar[] pwdDataset=new SVLVar[pwdTargetLigandObject.getNumberOfMembersInFile()+8];
            if(pwdDataset.length!=11)
            {
                return new SVLVar(new SVLVar(), new SVLVar("pwd set wrong size: " + pwdDataset.length + ".", true));
                //throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            }
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
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
    }
    else
    {
    SVLVar[] data=new SVLVar[pwdGroup.getNumberOfMembersInFile()];
    
        for(int index=0;index<pwdGroup.getNumberOfMembersInFile();index++)
        {
                    List<HObject> ligandList=pwdGroup.getMemberList();
                    H5Group pwdTargetLigandObject=(H5Group)ligandList.get(index);
            SVLVar[] pwdDataset=new SVLVar[pwdTargetLigandObject.getNumberOfMembersInFile()+8];
            if(pwdDataset.length!=11)
            {
                return new SVLVar(new SVLVar(), new SVLVar("pwd set wrong size: " + pwdDataset.length + ".", true));
                //throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            }
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
                }
                else if(member.getName().compareTo("Index B")==0)
                {
                    ligandIndex=(int[])member.read();
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
                        targetIndex[atomIndex]++;
                        ligandIndex[atomIndex]++;
                    }
                    pwdDataset[1]=new SVLVar(targetIndex);
                    pwdDataset[2]=new SVLVar(ligandIndex);
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
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
    }
}
    
private SVLVar retrieveAtomByAtomDecomposition(SVLVar var) throws SVLJavaException, IOException, Exception
{
                sessionErrBuffer=new PrintStream(new java.io.ByteArrayOutputStream(), true);
                sessionOutBuffer=new PrintStream(new java.io.ByteArrayOutputStream(), true);
        java.lang.System.setErr(sessionErrBuffer);
        java.lang.System.setOut(sessionOutBuffer);
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    String choice = var.peek(1).getTokn(3);
    String ligand="";
    if(var.peek(1).length()>3)
    {
        ligand = var.peek(1).getTokn(4);
        //if(true) return new SVLVar("retrieveAtomByAtomDecomposition "+ligand);
    }
//                    FileWriter pwdOutputFile=new FileWriter("./check.pwd");
//                    BufferedWriter bufferedPWDWriter=new BufferedWriter(pwdOutputFile);
//                        bufferedPWDWriter.append(target+" here\n");
//                        bufferedPWDWriter.flush();
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
                    String xPath="/DivCon/"+target+"/Pairwise Decomposition";
                    H5Group pwdGroup=(H5Group)findHDF5Object(h5File, xPath);
                    DivconType divcon=loadTarget(h5File, target);
    if(ligand.length()>0)
    {
                    ArrayList<Double> xsList=new ArrayList<>();
                    ArrayList<Double> ysList=new ArrayList<>();
                    ArrayList<Double> zsList=new ArrayList<>();
                    getTargetCoordinates(divcon, xsList, ysList, zsList);
                    double[] x=new double[xsList.size()];
                    double[] y=new double[ysList.size()];
                    double[] z=new double[zsList.size()];
                    for(int vIndex=0;vIndex<x.length;vIndex++)
                    {
                        x[vIndex]=xsList.get(vIndex);
                        y[vIndex]=ysList.get(vIndex);
                        z[vIndex]=zsList.get(vIndex);
                    }
//                    xPath="/DivCon/"+target+"/"+target;
//                    H5CompoundDS hTargetCollectionObject=(H5CompoundDS)findHDF5Object(h5File, xPath);
//                    hTargetCollectionObject.selectMember(0);
//                    String[] targetAtomSymbols=(String[])((Vector)hTargetCollectionObject.read()).elementAt(0);
//                    hTargetCollectionObject.selectMember(1);
//                    String[] targetAtomNames=(String[])((Vector)hTargetCollectionObject.read()).elementAt(1);
//                    hTargetCollectionObject.selectMember(2);
//                    String[] targetResidueNames=(String[])((Vector)hTargetCollectionObject.read()).elementAt(2);
//                    hTargetCollectionObject.selectMember(5);
//                    int[] targetSequence=(int[])((Vector)hTargetCollectionObject.read()).elementAt(5);
//                    hTargetCollectionObject.selectMember(7);
//                    double[] x=(double[])((Vector)hTargetCollectionObject.read()).elementAt(7);
//                    hTargetCollectionObject.selectMember(8);
//                    double[] y=(double[])((Vector)hTargetCollectionObject.read()).elementAt(8);
//                    hTargetCollectionObject.selectMember(9);
//                    double[] z=(double[])((Vector)hTargetCollectionObject.read()).elementAt(9);
//                    hTargetCollectionObject.selectMember(12);
//                    double[] epsilon=(double[])((Vector)hTargetCollectionObject.read()).elementAt(12);
//                    hTargetCollectionObject.selectMember(13);
//                    double[] welllDepth=(double[])((Vector)hTargetCollectionObject.read()).elementAt(13);
//                    hTargetCollectionObject.selectMember(14);
//                    int[] formalCharge=(int[])((Vector)hTargetCollectionObject.read()).elementAt(14);
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
//                    ListIterator li=hTargetCollectionObject.getMetadata().listIterator();
//                    while(li.hasNext())
//                    {
//                        Object obj=li.next();
//                        if(obj instanceof hdf.object.Attribute && ((hdf.object.Attribute)obj).getName().compareTo("Total Charge")==0)
//                        {
////                            ChargesType chargeType=new ChargesType();
////                            TotalChargeType totalChargeType=new TotalChargeType();
////                            chargeType.setTotalCharge(totalChargeType);
////                            int[] tc=(int[])((hdf.object.Attribute)obj).getValue();
//                        }
//                    }
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
            return new SVLVar(new SVLVar(), new SVLVar("", true));
            }
            SVLVar[] pwdDataset=new SVLVar[pwdTargetLigandObject.getNumberOfMembersInFile()+2];
            if(pwdDataset.length!=11)
            {
                return new SVLVar(new SVLVar(), new SVLVar("pwd set wrong size: " + pwdDataset.length + ".", true));
                //throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            }
            List pwdRow=pwdTargetLigandObject.getMemberList();
            pwdDataset[0]=new SVLVar(pwdTargetLigandObject.getName());
            int[] targetIndex=null;
            int[] ligandIndex=null;
            double[] pwdValues=null;
            H5CompoundDS pwdObject=null;
            for(int memberIndex=0;memberIndex<3;memberIndex++)
            {
                H5CompoundDS member=(H5CompoundDS)pwdRow.get(memberIndex);
                if(member.getName().compareTo("Atom A")==0)
                {
                    targetIndex=(int[])member.read();
                    pwdDataset[1]=new SVLVar((int[])member.read());
                }
                else if(member.getName().compareTo("Atom B")==0)
                {
                    ligandIndex=(int[])member.read();
                    pwdDataset[2]=new SVLVar((int[])member.read());
                }
                else
                {
                    pwdObject=(H5CompoundDS)member;
                    pwdValues=(double[])member.read();
                }
            }
//                    xPath="/DivCon/"+pwdObject.getName();
//                    H5CompoundDS hLigandCollectionObject=(H5CompoundDS)findHDF5Object(h5File, xPath);
//                    hLigandCollectionObject.selectMember(0);
//                    String[] ligandAtomSymbols=(String[])((Vector)hLigandCollectionObject.read()).elementAt(0);
//                    hLigandCollectionObject.selectMember(1);
//                    String[] ligandAtomNames=(String[])((Vector)hLigandCollectionObject.read()).elementAt(1);
//                    hLigandCollectionObject.selectMember(2);
//                    String[] ligandResidueNames=(String[])((Vector)hLigandCollectionObject.read()).elementAt(2);
//                    hLigandCollectionObject.selectMember(5);
//                    int[] ligandSequence=(int[])((Vector)hLigandCollectionObject.read()).elementAt(5);
//                    hLigandCollectionObject.selectMember(7);
//                    double[] xLigand=(double[])((Vector)hLigandCollectionObject.read()).elementAt(7);
//                    hLigandCollectionObject.selectMember(8);
//                    double[] yLigand=(double[])((Vector)hLigandCollectionObject.read()).elementAt(8);
//                    hLigandCollectionObject.selectMember(9);
//                    double[] zLigand=(double[])((Vector)hLigandCollectionObject.read()).elementAt(9);
//                    hLigandCollectionObject.selectMember(7);
//                    double[] epsilon=(double[])((Vector)hLigandCollectionObject.read()).elementAt(7);
//                    hLigandCollectionObject.selectMember(8);
//                    double[] welllDepth=(double[])((Vector)hLigandCollectionObject.read()).elementAt(8);
//                    hLigandCollectionObject.selectMember(9);
//                    int[] formalCharge=(int[])((Vector)hLigandCollectionObject.read()).elementAt(9);
//                    SVLJava.print(targetAtomSymbols.length+" lengths "+ligandAtomSymbols.length);
                    ArrayList<Double> EabList=new  ArrayList<Double>();
//                    ArrayList<Double> EabpList=new  ArrayList<Double>();
//                    ArrayList<Double> EabcList=new  ArrayList<Double>();
//                    ArrayList<Double> LennardJonesList=new  ArrayList<Double>();
//                    ArrayList<Double> dispersionList=new  ArrayList<Double>();
//                    ArrayList<Double> repulsionList=new  ArrayList<Double>();
//                    ArrayList<Double> electrostaticList=new  ArrayList<Double>();
                    ArrayList<Double> distanceList=new  ArrayList<Double>();
                    
//                        bufferedPWDWriter.append(targetIndex.length+" "+ligandIndex.length+"\n");
                    for(int atomIndex=0;atomIndex<targetIndex.length;atomIndex++)
                    {
                        int ligandOffsetIndex=ligandIndex[atomIndex];
//                        bufferedPWDWriter.append(targetIndex[atomIndex]+" "+ligandOffsetIndex+"\n");
                        ligandOffsetIndex-=x.length;
                        double distance=0.0;//Math.sqrt(Math.pow(x[targetIndex[atomIndex]]-xLigand[ligandOffsetIndex],2)+Math.pow(y[targetIndex[atomIndex]]-yLigand[ligandOffsetIndex],2)+Math.pow(z[targetIndex[atomIndex]]-zLigand[ligandOffsetIndex],2));
                        //if(distance<=4.0 || pwdValues[7*iw]!=0.0)
                        EabList.add(pwdValues[7*atomIndex]);
//                        EabpList.add(pwdValues[7*atomIndex+1]);
//                        EabcList.add(pwdValues[7*atomIndex+2]);
//                        LennardJonesList.add(pwdValues[7*atomIndex+3]);
//                        dispersionList.add(pwdValues[7*atomIndex+4]);
//                        repulsionList.add(pwdValues[7*atomIndex+5]);
//                        electrostaticList.add(pwdValues[7*atomIndex+6]);
                        distanceList.add(distance);
                    }
                    double[] Eab=new double[pwdDataset[1].length()];
//                    double[] Eabp=new double[pwdDataset[1].length()];
//                    double[] Eabc=new double[pwdDataset[1].length()];
//                    double[] LennardJones=new double[pwdDataset[1].length()];
//                    double[] dispersion=new double[pwdDataset[1].length()];
//                    double[] repulsion=new double[pwdDataset[1].length()];
//                    double[] electrostatic=new double[pwdDataset[1].length()];
                    double[] distance=new double[pwdDataset[1].length()];
                    for(int vIndex=0;vIndex<pwdDataset[1].length();vIndex++)
                    {
                    Eab[vIndex]=EabList.get(vIndex);
//                    Eabp[vIndex]=EabpList.get(vIndex);
//                    Eabc[vIndex]=EabcList.get(vIndex);
//                    LennardJones[vIndex]=LennardJonesList.get(vIndex);
//                    dispersion[vIndex]=dispersionList.get(vIndex);
//                    repulsion[vIndex]=repulsionList.get(vIndex);
//                    electrostatic[vIndex]=electrostaticList.get(vIndex);
                    distance[vIndex]=distanceList.get(vIndex);
                    }
                    pwdDataset[3]=new SVLVar(Eab);
//                    pwdDataset[4]=new SVLVar(Eabp);
//                    pwdDataset[5]=new SVLVar(Eabc);
//                    pwdDataset[6]=new SVLVar(LennardJones);
//                    pwdDataset[7]=new SVLVar(dispersion);
//                    pwdDataset[8]=new SVLVar(repulsion);
//                    pwdDataset[9]=new SVLVar(electrostatic);
                    pwdDataset[4]=new SVLVar(distance);
            SVLVar data=new SVLVar(new String[]{"name", "indexA", "indexB", "Eab", "distance"}, pwdDataset);
        h5File.close();
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
    }
    else
    {
        SVLVar[] data=new SVLVar[1];
        xPath="/DivCon/"+target+"/Pairwise Decomposition/"+choice;
        H5CompoundDS pwdData=(H5CompoundDS)findHDF5Object(h5File, xPath);
                    pwdData.selectMember(0);
                    int[] indexA=(int[])((Vector)pwdData.read()).elementAt(0);
                    pwdData.selectMember(1);
                    int[] indexB=(int[])((Vector)pwdData.read()).elementAt(1);
                    pwdData.selectMember(2);
                    double[] distance=(double[])((Vector)pwdData.read()).elementAt(2);
                    pwdData.selectMember(3);
                    double[] Eab=(double[])((Vector)pwdData.read()).elementAt(3);
            SVLVar[] pwdDataset=new SVLVar[5];
                    pwdDataset[0]=new SVLVar(pwdGroup.getName());
//                    double[] distance=new double[indexA.length];
                    for(int vIndex=0;vIndex<indexA.length;vIndex++)
                    {
//                        distance[vIndex]=Math.sqrt(Math.pow(x[indexA[vIndex]]-x[indexB[vIndex]],2)+Math.pow(y[indexA[vIndex]]-y[indexB[vIndex]],2)+Math.pow(z[indexA[vIndex]]-z[indexB[vIndex]],2));
                        indexA[vIndex]++;
                        indexB[vIndex]++;
                    }
                    pwdDataset[1]=new SVLVar(indexA);
                    pwdDataset[2]=new SVLVar(indexB);
                    pwdDataset[3]=new SVLVar(Eab);
                    pwdDataset[4]=new SVLVar(distance);
            data[0]=new SVLVar(new String[]{"name", "indexA", "indexB", "Eab", "distance"}, pwdDataset);
    
//        for(int index=0;index<pwdGroup.getNumberOfMembersInFile();index++)
//        {
////                    List<HObject> ligandList=pwdGroup.getMemberList();
////                    H5Group pwdTargetLigandObject=(H5Group)ligandList.get(index);
//            SVLVar[] pwdDataset=new SVLVar[pwdGroup.getNumberOfMembersInFile()+2];
//            if(pwdDataset.length!=4)
//              {
//                  return new SVLVar(new SVLVar(), new SVLVar("pwd set wrong size: " + pwdDataset.length + "."));
//                  //throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
//              }
//            List pwdRow=pwdGroup.getMemberList();
//            pwdDataset[0]=new SVLVar(pwdGroup.getName());
//            int[] targetIndex=null;
//            int[] ligandIndex=null;
//            double[] pwdValues=null;
//            H5CompoundDS pwdObject=null;
//            for(int memberIndex=0;memberIndex<2;memberIndex++)
//            {
//                H5CompoundDS member=(H5CompoundDS)pwdRow.get(memberIndex);
//                if(member.getName().compareTo("Atom A")==0)
//                {
//                    member.selectMember(0);
//                    targetIndex=(int[])((Vector)member.read()).elementAt(0);
//                    pwdDataset[1]=new SVLVar(targetIndex);
//                }
//                else if(member.getName().compareTo("Atom B")==0)
//                {
//                    member.selectMember(1);
//                    ligandIndex=(int[])((Vector)member.read()).elementAt(1);
//                    ligandIndex=(int[])member.read();
//                    pwdDataset[2]=new SVLVar((int[])member.read());
//                }
//                else
//                {
//                    member.selectMember(2);
//                    pwdObject=(H5CompoundDS)member;
//                    pwdValues=(double[])((Vector)member.read()).elementAt(2);
//                }
//            }
//                    /*xPath="/DivCon/"+pwdObject.getName();
//                    H5CompoundDS hLigandCollectionObject=(H5CompoundDS)findHDF5Object(h5File, xPath);
//                    hLigandCollectionObject.selectMember(0);
//                    String[] ligandAtomSymbols=(String[])((Vector)hLigandCollectionObject.read()).elementAt(0);
//                    hLigandCollectionObject.selectMember(1);
//                    String[] ligandAtomNames=(String[])((Vector)hLigandCollectionObject.read()).elementAt(1);
//                    hLigandCollectionObject.selectMember(2);
//                    String[] ligandResidueNames=(String[])((Vector)hLigandCollectionObject.read()).elementAt(2);
//                    hLigandCollectionObject.selectMember(5);
//                    int[] ligandSequence=(int[])((Vector)hLigandCollectionObject.read()).elementAt(5);
//                    hLigandCollectionObject.selectMember(7);
//                    double[] xLigand=(double[])((Vector)hLigandCollectionObject.read()).elementAt(7);
//                    hLigandCollectionObject.selectMember(8);
//                    double[] yLigand=(double[])((Vector)hLigandCollectionObject.read()).elementAt(8);
//                    hLigandCollectionObject.selectMember(9);
//                    double[] zLigand=(double[])((Vector)hLigandCollectionObject.read()).elementAt(9);*/
////                    hLigandCollectionObject.selectMember(7);
////                    double[] epsilon=(double[])((Vector)hLigandCollectionObject.read()).elementAt(7);
////                    hLigandCollectionObject.selectMember(8);
////                    double[] welllDepth=(double[])((Vector)hLigandCollectionObject.read()).elementAt(8);
////                    hLigandCollectionObject.selectMember(9);
////                    int[] formalCharge=(int[])((Vector)hLigandCollectionObject.read()).elementAt(9);
//                    //SVLJava.print(targetAtomSymbols.length+" lengths "+ligandAtomSymbols.length);
//                    ArrayList<Double> EabList=new  ArrayList<Double>();
////                    ArrayList<Double> EabpList=new  ArrayList<Double>();
////                    ArrayList<Double> EabcList=new  ArrayList<Double>();
////                    ArrayList<Double> LennardJonesList=new  ArrayList<Double>();
////                    ArrayList<Double> dispersionList=new  ArrayList<Double>();
////                    ArrayList<Double> repulsionList=new  ArrayList<Double>();
////                    ArrayList<Double> electrostaticList=new  ArrayList<Double>();
//                    ArrayList<Double> distanceList=new  ArrayList<Double>();
//                    
////                        bufferedPWDWriter.append(targetIndex.length+" "+ligandIndex.length+"\n");
//                    for(int atomIndex=0;atomIndex<targetIndex.length;atomIndex++)
//                    {
//                        int ligandOffsetIndex=ligandIndex[atomIndex];
////                        bufferedPWDWriter.append(targetIndex[atomIndex]+" "+ligandOffsetIndex+"\n");
//                        ligandOffsetIndex-=x.length;
//                        double distance=0.0;//Math.sqrt(Math.pow(x[targetIndex[atomIndex]]-xLigand[ligandOffsetIndex],2)+Math.pow(y[targetIndex[atomIndex]]-yLigand[ligandOffsetIndex],2)+Math.pow(z[targetIndex[atomIndex]]-zLigand[ligandOffsetIndex],2));
//                        //if(distance<=4.0 || pwdValues[7*iw]!=0.0)
//                        EabList.add(pwdValues[7*atomIndex]);
////                        EabpList.add(pwdValues[7*atomIndex+1]);
////                        EabcList.add(pwdValues[7*atomIndex+2]);
////                        LennardJonesList.add(pwdValues[7*atomIndex+3]);
////                        dispersionList.add(pwdValues[7*atomIndex+4]);
////                        repulsionList.add(pwdValues[7*atomIndex+5]);
////                        electrostaticList.add(pwdValues[7*atomIndex+6]);
//                        distanceList.add(distance);
//                    }
//                    double[] Eab=new double[pwdDataset[1].length()];
////                    double[] Eabp=new double[pwdDataset[1].length()];
////                    double[] Eabc=new double[pwdDataset[1].length()];
////                    double[] LennardJones=new double[pwdDataset[1].length()];
////                    double[] dispersion=new double[pwdDataset[1].length()];
////                    double[] repulsion=new double[pwdDataset[1].length()];
////                    double[] electrostatic=new double[pwdDataset[1].length()];
//                    double[] distance=new double[pwdDataset[1].length()];
//                    for(int vIndex=0;vIndex<pwdDataset[1].length();vIndex++)
//                    {
//                    Eab[vIndex]=EabList.get(vIndex);
////                    Eabp[vIndex]=EabpList.get(vIndex);
////                    Eabc[vIndex]=EabcList.get(vIndex);
////                    LennardJones[vIndex]=LennardJonesList.get(vIndex);
////                    dispersion[vIndex]=dispersionList.get(vIndex);
////                    repulsion[vIndex]=repulsionList.get(vIndex);
////                    electrostatic[vIndex]=electrostaticList.get(vIndex);
//                    distance[vIndex]=distanceList.get(vIndex);
//                    }
//                    pwdDataset[3]=new SVLVar(Eab);
////                    pwdDataset[4]=new SVLVar(Eabp);
////                    pwdDataset[5]=new SVLVar(Eabc);
////                    pwdDataset[6]=new SVLVar(LennardJones);
////                    pwdDataset[7]=new SVLVar(dispersion);
////                    pwdDataset[8]=new SVLVar(repulsion);
////                    pwdDataset[9]=new SVLVar(electrostatic);
//                    pwdDataset[4]=new SVLVar(distance);
//            data[index]=new SVLVar(new String[]{"name", "indexA", "indexB", "Eab", "distance"}, pwdDataset);
//        }
//             bufferedPWDWriter.close();       
        h5File.close();
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
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
                        if(obj instanceof hdf.object.Attribute && ((hdf.object.Attribute)obj).getName().compareTo("Total Charge")==0)
                        {
//                            ChargesType chargeType=new ChargesType();
//                            TotalChargeType totalChargeType=new TotalChargeType();
//                            chargeType.setTotalCharge(totalChargeType);
//                            int[] tc=(int[])((hdf.object.Attribute)obj).getValue();
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
            return new SVLVar(new SVLVar(), new SVLVar("", true));
            }
            SVLVar[] pwdDataset=new SVLVar[pwdTargetLigandObject.getNumberOfMembersInFile()+4];
            if(pwdDataset.length!=7)
            {
                return new SVLVar(new SVLVar(), new SVLVar("pwd set wrong size: " + pwdDataset.length + ".", true));
                //throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            }
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
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
    }
    else
    {
        SVLVar[] data=new SVLVar[pwdGroup.getNumberOfMembersInFile()];
    
        for(int index=0;index<pwdGroup.getNumberOfMembersInFile();index++)
        {
                    List<HObject> ligandList=pwdGroup.getMemberList();
                    H5Group pwdTargetLigandObject=(H5Group)ligandList.get(index);
            SVLVar[] pwdDataset=new SVLVar[pwdTargetLigandObject.getNumberOfMembersInFile()+4];
            if(pwdDataset.length!=7)
            {
                return new SVLVar(new SVLVar(), new SVLVar("pwd set wrong size: " + pwdDataset.length + ".", true));
                //throw new SVLJavaException("pwd set wrong size: " + pwdDataset.length + ".");
            }
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
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
    }
}
    
private SVLVar retrieveNMRScore(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    String ligand = "";
    if(var.peek(1).length()>2)
    {
        ligand = var.peek(1).getTokn(3);
        //if(true) return new SVLVar("retrieveAtomByAtomPWD "+ligand);
    }
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
        String xPath="/DivCon/"+target+"/NMR Score Class Number Average";
        HObject ho=findHDF5Object(h5File, xPath);
        H5CompoundDS hObject=(H5CompoundDS)ho;
        if(hObject==null)
        {
            return new SVLVar(new SVLVar(""), new SVLVar("'NMR Score Class Number Average' does not exist for: '" + target + "'.", true));
            //throw new SVLJavaException("'NMR Score Class Number Average' does not exist for: '" + target + "'.");
        }
        hObject.init();
        hObject.clear();
        hObject.setMemberSelection(true);
        Vector o=(Vector)hObject.getData();
//        String[] dimNames=hObject.getDimNames();
        SVLVar[] data=new SVLVar[hObject.getMemberCount()];
        if(ligand.isEmpty())
        {
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
        }
        else
        {
            int rowIndex=0;
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
                    boolean hit=false;
                    while(!hit && rowIndex<rowData.length)
                    {
                        if(rowData[rowIndex].compareTo(ligand)==0)
                        {
                            hit=true;
                        }
                        else
                        {
                            rowIndex++;
                        }
                    }
                    if(!hit)new SVLVar(new SVLVar(""), new SVLVar(ligand+" not in 'NMR Score Class Number Average' table for '" + target + "'.", true));
                    
                    data[index]=new SVLVar(new String[]{rowData[rowIndex]});
    //                for(int columnIndex=0;columnIndex<rowData.length;columnIndex++)
    //                {
    //                    //scores[columnIndex]=(double)rowData[columnIndex];
    //                }
                }
                else
                {
                    //double[] scores = new double[11];
                    double[] rowData=(double[])o.elementAt(index);
                    data[index]=new SVLVar(new double[]{rowData[rowIndex]});
    //                for(int columnIndex=0;columnIndex<rowData.length;columnIndex++)
    //                {
    //                    scores[columnIndex]=(double)rowData[columnIndex];
    //                }
                }
            }
        }
            h5File.close();
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
}
    
private SVLVar retrieveChemicalShifts(SVLVar var) throws SVLJavaException, IOException, Exception
{
    sessionErrBuffer=new PrintStream(new java.io.ByteArrayOutputStream(), true);
    sessionOutBuffer=new PrintStream(new java.io.ByteArrayOutputStream(), true);
    java.lang.System.setErr(sessionErrBuffer);
    java.lang.System.setOut(sessionOutBuffer);
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    String ligand="";
    if(var.peek(1).length()>2)
    {
        ligand = var.peek(1).getTokn(3);
        //if(true) return new SVLVar("retrieveAtomByAtomPWD "+ligand);
    }
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
        String xPath="/DivCon/"+target+"/NMR Matrices";
        HObject ho=findHDF5Object(h5File, xPath);
        if(ho==null) return new SVLVar(new SVLVar(), new SVLVar("", true));
    if(ligand.length()>0)
    {
        xPath="/DivCon/"+target+"/NMR Matrices/"+ligand;
        ho=findHDF5Object(h5File, xPath);
        H5CompoundDS hObject=(H5CompoundDS)ho;
        if(hObject==null)
        {
            return new SVLVar(new SVLVar(), new SVLVar("NMR Score does not exist for: '" + target + "'.", true));
            //throw new SVLJavaException("NMR Score does not exist for: '" + target + "'.");
        }
        hObject.init();
        hObject.clear();
        hObject.setMemberSelection(true);
        Vector o=(Vector)hObject.getData();
        String[] dimNames=hObject.getDimNames();
        SVLVar[] chemicalShiftData=new SVLVar[hObject.getMemberCount()];
        for(int index=0;index<hObject.getMemberCount();index++)
        {
            hObject.selectMember(index);
            if(index<=1)
            {
                int[] indexData=(int[])o.elementAt(index);
                if(index==0)
                for(int indexIndex=0;indexIndex<indexData.length;indexIndex++)
                {
                    indexData[indexIndex]+=1;
                }
                chemicalShiftData[index]=new SVLVar(indexData);
            }
            else
            {
                //double[] scores = new double[11];
                double[] rowData=(double[])o.elementAt(index);
                chemicalShiftData[index]=new SVLVar(rowData);
//                for(int columnIndex=0;columnIndex<rowData.length;columnIndex++)
//                {
//                    scores[columnIndex]=(double)rowData[columnIndex];
//                }
            }
        }
        
           SVLVar data=new SVLVar(new String[]{"Index", "Role", "Bound Shift", "Unbound Shift", "Exp Bound Shift", "Exp Unbound Shift"}, chemicalShiftData);
            h5File.close();
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
    }
    else
    {
    SVLVar[] data=new SVLVar[((H5Group)ho).getNumberOfMembersInFile()];
//        List pwdMembers=hObject.getMemberList();
//        for(int index=0;index<hObject.getNumberOfMembersInFile();index++)
//        {
//            H5Group pwdGroup=(H5Group)pwdMembers.get(index);
//            SVLVar[] chemicalShiftsDataset=new SVLVar[pwdGroup.getNumberOfMembersInFile()+1];
//            if(chemicalShiftsDataset.length!=5) throw new SVLJavaException("pwd set wrong size: " + chemicalShiftsDataset.length + ".");
//            List pwdRow=pwdGroup.getMemberList();
//            chemicalShiftsDataset[0]=new SVLVar(pwdGroup.getName());
//            for(int memberIndex=0;memberIndex<4;memberIndex++)
//            {
//                H5ScalarDS member=(H5ScalarDS)pwdRow.get(memberIndex);
//                if(member.getName().compareTo("Bound Shift")==0)
//                {
//                    chemicalShiftsDataset[1]=new SVLVar((int[])member.read());
//                }
//                else if(member.getName().compareTo("Unbound Shift")==0)
//                {
//                    chemicalShiftsDataset[2]=new SVLVar((int[])member.read());
//                }
//                else if(member.getName().compareTo("Exp Bound Shift")==0)
//                {
//                    chemicalShiftsDataset[3]=new SVLVar(member.readBytes());
//                }
//                else
//                {
//                    chemicalShiftsDataset[4]=new SVLVar((double[])member.read());
//                }
//            }
////            pwdDataset[1]=new SVLVar((int[])sequenceA.read());
////            H5ScalarDS sequenceB=(H5ScalarDS)pwdRow.get(1);
////            pwdDataset[2]=new SVLVar((int[])sequenceB.read());
////            H5ScalarDS sequenceValues=(H5ScalarDS)pwdRow.get(2);
////            pwdDataset[3]=new SVLVar((double[])sequenceValues.read());
////            H5ScalarDS sequenceLabels=(H5ScalarDS)pwdRow.get(3);
////            pwdDataset[4]=new SVLVar(sequenceLabels.readBytes());
//            data[index]=new SVLVar(new String[]{"name", "Bound Shift", "Unbound Shift", "Exp Bound Shift", "Exp Unbound Shift"}, chemicalShiftsDataset);
//        }
            h5File.close();
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
    }
}
    
private SVLVar retrieveDensities(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
    String xPath="/DivCon/"+target+"/Densities";
    H5Group nmrGroup=(H5Group)findHDF5Object(h5File, xPath);
    SVLVar[] averages=new SVLVar[3];
    if(nmrGroup!=null)
    {
        List nmrRow=nmrGroup.getMemberList();
        for(int memberIndex=0;memberIndex<1;memberIndex++)
        {
            H5ScalarDS member=(H5ScalarDS)nmrRow.get(memberIndex);
            double[] densities=(double[])member.read();
            averages[0]=new SVLVar(Math.sqrt(densities.length));
            averages[1]=new SVLVar(Math.sqrt(densities.length));
            averages[2]=new SVLVar(densities);
        }
    }
    else
    {
        averages[0]=new SVLVar(0);
        averages[1]=new SVLVar(0);
        averages[2]=new SVLVar();
    }
    SVLVar data=new SVLVar(new String[]{"rows", "cols", "densities"}, averages);
            h5File.close();
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
}
    
private SVLVar retrieveEigenVectors(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
    String xPath="/DivCon/"+target+"/Eigenvectors";
    H5Group nmrGroup=(H5Group)findHDF5Object(h5File, xPath);
    SVLVar[] averages=new SVLVar[3];
    if(nmrGroup!=null)
    {
        List nmrRow=nmrGroup.getMemberList();
        for(int memberIndex=0;memberIndex<1;memberIndex++)
        {
            H5ScalarDS member=(H5ScalarDS)nmrRow.get(memberIndex);
            double[] densities=(double[])member.read();
            averages[0]=new SVLVar(Math.sqrt(densities.length));
            averages[1]=new SVLVar(Math.sqrt(densities.length));
            averages[2]=new SVLVar(densities);
        }
    }
    else
    {
        averages[0]=new SVLVar(0);
        averages[1]=new SVLVar(0);
        averages[2]=new SVLVar();
    }
    SVLVar data=new SVLVar(new String[]{"rows", "cols", "eigenVectors"}, averages);
            h5File.close();
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
}
    
private SVLVar retrieveEnergyLevels(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
//   java.lang.System.out.println("h5 opened ");
    String xPath="/DivCon/"+target+"/Energy Levels";
    H5Group nmrGroup=(H5Group)findHDF5Object(h5File, xPath);
    SVLVar[] averages=new SVLVar[1];
    if(nmrGroup!=null)
    {
        List nmrRow=nmrGroup.getMemberList();
        for(int memberIndex=0;memberIndex<1;memberIndex++)
        {
            H5ScalarDS member=(H5ScalarDS)nmrRow.get(memberIndex);
            averages[0]=new SVLVar((double[])member.read());
        }
    }
    else
    {
        averages[0]=new SVLVar(new String[]{filename, target});
    }
    SVLVar data=new SVLVar(new String[]{"energyLevels"}, averages);
    h5File.close();
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
}
    
    private SVLVar setNMRAtomSelection(SVLVar var) throws SVLJavaException, IOException, Exception
    {
        String filename = var.peek(1).getTokn(1);
        String target = var.peek(1).getTokn(2);
        H5File h5File=new H5File(filename, H5File.WRITE);
        h5File.open();
            String xPath="/DivCon/"+target+"/Selections";
            HObject ho=findHDF5Object(h5File, xPath);
            H5Group hObject=(H5Group)ho;
            if(hObject!=null){
                H5Group selectionGroup=(H5Group)h5File.createGroup("Selection", hObject);
                long[] dims=new long[1];
                dims[0]=var.peek(1).getInt(3);
                H5Datatype h5Datatype=new H5Datatype(H5Datatype.CLASS_INTEGER);
                H5ScalarDS member=(H5ScalarDS)h5File.createScalarDS("NMR Selection", selectionGroup, h5Datatype, dims, null, null, 0, 0, var.peek(1).peek(4).getInts());
                long[] attrDims = { 1 };
                 Attribute attr = new Attribute("Role", new H5Datatype(Datatype.CLASS_STRING, 5, -1, -1), attrDims);
                 attr.setValue("Atoms");
 
                h5File.writeAttribute(member, attr, true);
            }
        h5File.close();
        return new SVLVar(new SVLVar(), new SVLVar("",true));
    }
    
    protected int indexW(int i, int j)
    {
        int indw,ki,kj;
        ki=Math.max(i,j);
        kj=Math.min(i,j);
        indw=(ki*(ki-1))/2+kj;
        return indw;
    }
    
private SVLVar retrieveHamiltonian(SVLVar var) throws SVLJavaException, IOException, Exception
{
        String target = var.peek(1).getTokn(1);
    String hamiltonianParameterFile = var.peek(1).getTokn(2);
    String qbHome="/home/roger/testing-grounds/qmechanic";//java.lang.System.getProperty("QBHOME");
//    if(true)return new SVLVar(qbHome);
    File f=new File(qbHome+File.separator+"data"+File.separator+"Hamiltonian"+File.separator+hamiltonianParameterFile+".xml");
    JAXBContext jc = JAXBContext.newInstance("com.quantumbioinc.xml");
    Unmarshaller um=jc.createUnmarshaller();
    Hamiltonian hamiltonian=(Hamiltonian)um.unmarshal(f);
    ListIterator<Object> listIterator = hamiltonian.getHamiltonianParams().getElementAndPair().listIterator();
    int counter=0;
    SVLVar[] aCore=new SVLVar[83];
    SVLVar[] uss=new SVLVar[83];
    SVLVar[] upp=new SVLVar[83];
    SVLVar[] udd=new SVLVar[83];
    SVLVar[] zetas=new SVLVar[83];
    SVLVar[] zetap=new SVLVar[83];
    SVLVar[] zetad=new SVLVar[83];
    SVLVar[] zetasn=new SVLVar[83];
    SVLVar[] zetapn=new SVLVar[83];
    SVLVar[] zetadn=new SVLVar[83];
    while(counter<83 && listIterator.hasNext())
    {
        Object obj=listIterator.next();
        if(obj instanceof Element)
        {
            Element e=(Element)obj;
            if(e.getAcore()!=null){aCore[e.getNumber()-1]=new SVLVar(e.getAcore().getValue());}else{aCore[e.getNumber()-1]=new SVLVar();};
            if(e.getUss()!=null){uss[e.getNumber()-1]=new SVLVar(e.getUss().getValue());}else{uss[e.getNumber()-1]=new SVLVar();};
            if(e.getUpp()!=null){upp[e.getNumber()-1]=new SVLVar(e.getUpp().getValue());}else{upp[e.getNumber()-1]=new SVLVar();};
            if(e.getUdd()!=null){udd[e.getNumber()-1]=new SVLVar(e.getUdd().getValue());}else{udd[e.getNumber()-1]=new SVLVar();};
            if(e.getZetas()!=null){zetas[e.getNumber()-1]=new SVLVar(e.getZetas().getValue());}else{zetas[e.getNumber()-1]=new SVLVar();};
            if(e.getZetap()!=null){zetap[e.getNumber()-1]=new SVLVar(e.getZetap().getValue());}else{zetap[e.getNumber()-1]=new SVLVar();};
            if(e.getZetad()!=null){zetad[e.getNumber()-1]=new SVLVar(e.getZetad().getValue());}else{zetad[e.getNumber()-1]=new SVLVar();};
            if(e.getZetasn()!=null){zetasn[e.getNumber()-1]=new SVLVar(e.getZetasn().getValue());}else{zetasn[e.getNumber()-1]=new SVLVar();};
            if(e.getZetapn()!=null){zetapn[e.getNumber()-1]=new SVLVar(e.getZetapn().getValue());}else{zetapn[e.getNumber()-1]=new SVLVar();};
            if(e.getZetadn()!=null){zetadn[e.getNumber()-1]=new SVLVar(e.getZetadn().getValue());}else{zetadn[e.getNumber()-1]=new SVLVar();};
            counter++;
        }
    }
    for(int index=0;index<83;index++)
    {
        if(aCore[index]==null)aCore[index]=new SVLVar();
        if(uss[index]==null)uss[index]=new SVLVar();
        if(upp[index]==null)upp[index]=new SVLVar();
        if(udd[index]==null)udd[index]=new SVLVar();
        if(zetas[index]==null)zetas[index]=new SVLVar();
        if(zetap[index]==null)zetap[index]=new SVLVar();
        if(zetad[index]==null)zetad[index]=new SVLVar();
        if(zetasn[index]==null)zetasn[index]=new SVLVar();
        if(zetapn[index]==null)zetapn[index]=new SVLVar();
        if(zetadn[index]==null)zetadn[index]=new SVLVar();
    }
    SVLVar[] hamiltonianParameters=new SVLVar[9];
    hamiltonianParameters[0]=new SVLVar(uss);
    hamiltonianParameters[1]=new SVLVar(upp);
    hamiltonianParameters[2]=new SVLVar(udd);
    hamiltonianParameters[3]=new SVLVar(zetas);
    hamiltonianParameters[4]=new SVLVar(zetap);
    hamiltonianParameters[5]=new SVLVar(zetad);
    hamiltonianParameters[6]=new SVLVar(zetasn);
    hamiltonianParameters[7]=new SVLVar(zetapn);
    hamiltonianParameters[8]=new SVLVar(zetadn);

    SVLVar data=new SVLVar(new String[]{"Uss", "Upp", "Udd", "Zetas", "Zetap", "Zetad", "Zetasn", "Zetapn", "Zetadn"}, hamiltonianParameters);
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
}

private SVLVar retrieveDefaultProgramOptions(SVLVar var) throws SVLJavaException, IOException, Exception
{
            SVLVar[] defaultProgramOptions=new SVLVar[14];
            defaultProgramOptions[0]=new SVLVar("pm6", true);
            defaultProgramOptions[1]=new SVLVar("off", true);
            defaultProgramOptions[2]=new SVLVar("off", true);
            defaultProgramOptions[3]=new SVLVar("off", true);
            defaultProgramOptions[4]=new SVLVar("off", true);
            defaultProgramOptions[5]=new SVLVar(1.0e-8);
            defaultProgramOptions[6]=new SVLVar("off", true);
            defaultProgramOptions[7]=new SVLVar("off", true);
            defaultProgramOptions[8]=new SVLVar(new double[]{});
            defaultProgramOptions[9]=new SVLVar(new double[]{});
            defaultProgramOptions[10]=new SVLVar("off", true);
            defaultProgramOptions[11]=new SVLVar("off", true);
            defaultProgramOptions[12]=new SVLVar("off", true);
            defaultProgramOptions[13]=new SVLVar(1);
    return new SVLVar( new SVLVar(new String[]{"hamiltonian", "gradient", "opt", "freq", "decompose", "ecrit", "pwd",
                                   "perception", "selection", "region", "test", "restart", "standard", "np"}, defaultProgramOptions), new SVLVar("",true));
}
//[ '3FVA.pdb',
//[ ['3FVA.pdb.A','3FVA.pdb.W'], ['3FVA.pdb','3FVA.pdb'], ['',''], [6,2] ], 
//[ ['ASN','ASN','GLN','ASN','THR','PHE','HOH','HOH'], [1,2,3,4,5,6,7,8], "        ", ['amino','amino','amino','amino','amino','amino','none','none'], [16,14,17,14,14,21,3,3] ], 
//
//[ ['N','C','C','O','C','C','O','N','H','H','H','H','H','H','H','H','N','C','C','O','C','C','O','N','H','H','H','H','H','H','N','C','C','O','C','C','C','O','N','H','H','H','H','H','H','H','H','N','C','C','O','C','C','O','N','H','H','H','H','H','H','N','C','C','O','C','O','C','H','H','H','H','H','H','H','N','C','C','O','C','C','C','C','C','C','C','O','H','H','H','H','H','H','H','H','H','O','H','H','O','H','H'],
//  [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
//  ['sp3','sp3','sp2','sp2','sp3','sp2','sp2','sp2','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp2','sp3','sp2','sp2','sp3','sp2','sp2','sp2','sp3','sp3','sp3','sp3','sp3','sp3','sp2','sp3','sp2','sp2','sp3','sp3','sp2','sp2','sp2','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp2','sp3','sp2','sp2','sp3','sp2','sp2','sp2','sp3','sp3','sp3','sp3','sp3','sp3','sp2','sp3','sp2','sp2','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp2','sp3','sp2','sp2','sp3','sp2','sp2','sp2','sp2','sp2','sp2','sp2','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3','sp3'],
//  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
//  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
//  [ [2,9,10,11], [1,3,5,12], [2,4,17], 3, [2,6,13,14], [5,7,8], 6, [6,15,16], 1, 1, 1, 2, 5, 5, 8, 8, [3,18,25], [17,19,21,26], [18,20,31], 19, [18,22,27,28], [21,23,24], 22, [22,29,30], 17, 18, 21, 21, 24, 24, [19,32,40], [31,33,35,41], [32,34,48], 33, [32,36,42,43], [35,37,44,45], [36,38,39], 37, [37,46,47], 31, 32, 35, 35, 36, 36, 39, 39, [33,49,56], [48,50,52,57], [49,51,62], 50, [49,53,58,59], [52,54,55], 53, [53,60,61], 48, 49, 52, 52, 55, 55, [50,63,69], [62,64,66,70], [63,65,76], 64, [63,67,68,71], [66,72], [66,73,74,75], 62, 63, 66, 67, 68, 68, 68, [64,77,88], [76,78,80,89], [77,79,87], 78, [77,81,90,91], [80,82,83], [81,84,92], [81,85,93], [82,86,94], [83,86,95], [84,85,96], 78, 76, 77, 80, 80, 82, 83, 84, 85, 86, [98,99], 97, 97, [101,102], 100, 100 ],
//  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
//  ['N','CA','C','O','CB','CG','OD1','ND2','H1','H2','H3','HA','HB2','HB3','HD21','HD22','N','CA','C','O','CB','CG','OD1','ND2','H','HA','HB2','HB3','HD21','HD22','N','CA','C','O','CB','CG','CD','OE1','NE2','H','HA','HB2','HB3','HG2','HG3','HE21','HE22','N','CA','C','O','CB','CG','OD1','ND2','H','HA','HB2','HB3','HD21','HD22','N','CA','C','O','CB','OG1','CG2','H','HA','HB','HG1','HG21','HG22','HG23','N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ','OXT','H','HA','HB2','HB3','HD1','HD2','HE1','HE2','HZ','O','H1','H2','O','H1','H2'],
//  [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
//  [-6.702,-7.086,-6.105,-6.086,-8.507,-9.595,-9.322,-10.845,-7.491,-6.256,-6.113,-7.06,-8.6,-8.648,-11.502,-10.997,-5.305,-4.238,-4.339,-4.46,-2.891,-2.771,-2.94,-2.553,-5.363,-4.317,-2.797,-2.179,-2.478,-2.488,-4.312,-4.492,-3.406,-3.213,-5.873,-7.049,-8.384,-8.538,-9.375,-4.19,-4.427,-5.962,-5.937,-7.077,-6.928,-10.15,-9.238,-2.688,-1.748,-2.009,-2.262,-0.299,0.697,0.899,1.318,-2.727,-1.852,-0.126,-0.172,1.895,1.144,-1.942,-1.937,-0.871,-0.674,-3.322,-3.203,-3.844,-1.899,-1.718,-3.961,-3.928,-4.699,-3.94,-3.232,-0.187,0.767,0.043,-1.185,1.92,2.818,2.559,3.889,3.379,4.703,4.446,0.644,-0.256,1.156,1.55,2.454,1.84,4.063,3.205,5.421,4.999,-4.963,-4.55,-4.857,3.229,4.206,2.873],
//  [-0.246,-0.456,0.245,1.471,0.053,-0.823,-1.82,-0.455,-0.083,-1.029,0.509,-1.415,0.946,0.075,-0.909,0.237,-0.56,-0.059,-0.686,-1.902,-0.381,0.218,1.432,-0.63,-1.418,0.914,-1.343,-0.018,-0.334,-1.473,0.157,-0.3,0.301,1.513,0.128,-0.532,0.084,1.308,-0.766,1.005,-1.277,1.089,-0.111,-1.474,-0.424,-0.465,-1.614,-0.55,-0.074,-0.747,-1.949,-0.363,0.175,1.379,-0.718,-1.406,0.894,0.057,-1.322,-0.461,-1.555,0.035,-0.543,0.104,1.318,-0.437,-0.912,0.981,0.894,-1.495,-0.989,-0.862,1.011,1.305,1.55,-0.71,-0.2,0.142,-0.034,-1.196,-1.298,-2.231,-0.428,-2.309,-0.496,-1.433,0.627,-1.567,0.63,-2.076,-0.906,-2.814,0.209,-2.941,0.088,-1.488,1.946,2.566,2.342,1.122,1.113,1.272],
//  [-0.738,-2.167,-3.087,-3.179,-2.45,-1.836,-1.146,-2.097,-0.215,-0.408,-0.673,-2.368,-2.082,-3.409,-1.779,-2.584,-3.78,-4.637,-6.017,-6.14,-3.992,-2.611,-2.438,-1.607,-3.769,-4.728,-3.914,-4.542,-0.803,-1.763,-7.046,-8.432,-9.31,-9.314,-8.959,-8.226,-8.604,-8.617,-8.883,-6.972,-8.472,-8.857,-9.897,-8.454,-7.27,-9.102,-8.843,-10.038,-11.048,-12.39,-12.45,-10.629,-11.625,-11.723,-12.38,-9.967,-11.16,-9.772,-10.563,-12.963,-12.287,-13.461,-14.799,-15.676,-15.627,-15.495,-16.838,-15.511,-13.442,-14.724,-15.018,-17.218,-15.947,-14.612,-15.984,-16.472,-17.457,-18.765,-18.875,-17.666,-16.478,-15.482,-16.318,-14.355,-15.202,-14.213,-19.73,-16.463,-17.111,-17.839,-18.422,-15.571,-16.973,-13.695,-15.114,-13.467,-0.298,-0.935,0.592,-20.636,-20.715,-21.537] ] ]


private SVLVar retrieveLigandSelection(SVLVar var) throws SVLJavaException, IOException, Exception
{
    String filename = var.peek(1).getTokn(1);
    String target = var.peek(1).getTokn(2);
    H5File h5File=new H5File(filename, H5File.READ);
    h5File.open();
//   java.lang.System.out.println("h5 opened ");
    String xPath="/DivCon/"+target+"/Selections";
    H5Group selectionsGroup=(H5Group)findHDF5Object(h5File, xPath);
    SVLVar[] selections=new SVLVar[1];
    if(selectionsGroup!=null)
    {
        List nmrRow=selectionsGroup.getMemberList();
        for(int memberIndex=0;memberIndex<nmrRow.size();memberIndex++)
        {
            HObject hObject=(HObject)nmrRow.get(memberIndex);
            switch(hObject.getName())
            {
//                case "Compliment":
////            if(true) return new SVLVar(new SVLVar("here2", true), new SVLVar(""+((H5ScalarDS)hObject).getFullName(),true));
//                    try
//                    {
//                        selections[0]=new SVLVar((int[])((H5ScalarDS)hObject).read());
//                    }
//                    catch(HDF5Exception ex)
//                    {
//                        selections[0]=new SVLVar();
//                    }
//                     break;
                case "Wildtype":
//            if(true) return new SVLVar(new SVLVar("here", true), new SVLVar(""+((H5ScalarDS)hObject).getFullName(),true));
                    try
                    {
                        int[] indexData=(int[])((H5ScalarDS)hObject).read();
                        for(int indexIndex=0;indexIndex<indexData.length;indexIndex++)
                        {
                            indexData[indexIndex]+=1;
                        }
                        selections[0]=new SVLVar(indexData);
                    }
                    catch(HDF5Exception ex)
                    {
                        selections[0]=new SVLVar();
                    }
                    break;
            }
        }
    }
    else
    {
        return new SVLVar(new SVLVar("",true), new SVLVar("",true));
    }
    SVLVar data=new SVLVar(new String[]{"ligandSelection"}, selections);
    h5File.close();
    return new SVLVar(new SVLVar(data), new SVLVar("",true));
}
    private SVLVar retrieveTopologyClassNumbers(SVLVar var) throws SVLJavaException, IOException, Exception
    {
        String target = var.peek(1).getTokn(1);
                sessionErrBuffer=new PrintStream(new java.io.ByteArrayOutputStream(), true);
                sessionOutBuffer=new PrintStream(new java.io.ByteArrayOutputStream(), true);
        java.lang.System.setErr(sessionErrBuffer);
        java.lang.System.setOut(sessionOutBuffer);
        long[] dims = {1};
                    ObjectFactory objectFactory=new ObjectFactory();
                    DivconType divcon=new DivconType();
                    divcon.setVersion("1.0");
                    SVLVar svlMol=var.peek(1).peek(2);
                    String[] chainNames=svlMol.peek(2).peek(1).getTokns();
                    int[] residueCounts=svlMol.peek(2).peek(4).getInts();
                    String[] residueNames=svlMol.peek(3).peek(1).getTokns();
                    int[] sequences=svlMol.peek(3).peek(2).getInts();
                    int[] atomCounts=svlMol.peek(3).peek(5).getInts();
                    String[] symbols=svlMol.peek(4).peek(1).getTokns();
                    String[] hybridizations=svlMol.peek(4).peek(3).getTokns();
//                    int[] atomCounts=svlMol.peek(4).peek(5).getInts();
                    String[] atomNames=svlMol.peek(4).peek(8).getTokns();
                    int[] formalCharges=svlMol.peek(4).peek(2).getInts();
                    double[] x=svlMol.peek(4).peek(10).getReals();
                    double[] y=svlMol.peek(4).peek(11).getReals();
                    double[] z=svlMol.peek(4).peek(12).getReals();
                    Cml cml=objectFactory.createCml();
                    int chainIndex=0;
                    int residueIndex=0;
                    int chainOffset=0;
                    int residueOffset=0;
                    JAXBElement<Molecule> molecule=objectFactory.createMolecule(new Molecule());
                    JAXBElement<Molecule> submolecule=objectFactory.createMolecule(new Molecule());
                    JAXBElement<AtomArray> atomArray=objectFactory.createAtomArray(new AtomArray());
                    int totalCharge=0;
                    for(int atomIndex=0;atomIndex<atomNames.length;atomIndex++)
                    {
                        JAXBElement<Atom> atom=objectFactory.createAtom(new Atom());
                        atom.getValue().setElementType(symbols[atomIndex]);
                        atom.getValue().setX3(new Double(x[atomIndex]));
                        atom.getValue().setY3(new Double(y[atomIndex]));
                        atom.getValue().setZ3(new Double(z[atomIndex]));
                        atom.getValue().setFormalCharge(new BigInteger(""+formalCharges[atomIndex]));
                        totalCharge+=formalCharges[atomIndex];
                        JAXBElement<AtomType> atomType=objectFactory.createAtomType(new AtomType());
                        atomType.getValue().setName(atomNames[atomIndex]);
                        atom.getValue().getAnyCmlOrAnyOrAny().add(atomType);
                        JAXBElement<Scalar> chainID=objectFactory.createScalar(new Scalar());
                        chainID.getValue().setTitle("chainID");
                        chainID.getValue().setValue(chainNames[chainIndex].substring(chainNames[chainIndex].lastIndexOf('.')+1, chainNames[chainIndex].length()));
                        atom.getValue().getAnyCmlOrAnyOrAny().add(chainID);
                        JAXBElement<Scalar> hybridization=objectFactory.createScalar(new Scalar());
                        hybridization.getValue().setTitle("hybridization");
                        hybridization.getValue().setValue(hybridizations[atomIndex]);
                        atom.getValue().getAnyCmlOrAnyOrAny().add(hybridization);
                        atomArray.getValue().getAnyCmlOrAnyOrAny().add(atom);
                        if(atomIndex>=atomCounts[residueIndex]+residueOffset-1)
                        {
                            residueOffset+=atomCounts[residueIndex];
                            submolecule.getValue().setTitle(residueNames[residueIndex]);
                            submolecule.getValue().getAnyCmlOrAnyOrAny().add(atomArray);
                            JAXBElement<Scalar> sequence=objectFactory.createScalar(new Scalar());
                            sequence.getValue().setTitle("sequence");
                            sequence.getValue().setValue(""+sequences[residueIndex]);
                            submolecule.getValue().getAnyCmlOrAnyOrAny().add(sequence);
                            molecule.getValue().getAnyCmlOrAnyOrAny().add(submolecule);
                            residueIndex++;
                            if(residueIndex<=atomCounts.length)
                            {
                                atomArray=objectFactory.createAtomArray(new AtomArray());
                                submolecule=objectFactory.createMolecule(new Molecule());
                            }
                            if(residueIndex>=residueCounts[chainIndex]+chainOffset)
                            {
                                chainOffset+=residueCounts[chainIndex];
                                chainIndex++;
                                molecule.getValue().setTitle(svlMol.peek(1).getTokn(1));
                                cml.getAnyCmlOrAnyOrAny().add(molecule);
                                if(chainIndex<=residueCounts.length)molecule=objectFactory.createMolecule(new Molecule());
                            }
                        }
                    }
                    HamiltonianType hamiltonianType=objectFactory.createHamiltonianType();
                    hamiltonianType.setParameters("pm6");
                    divcon.setHamiltonian(hamiltonianType);
                    TotalChargeType totalChargeType=objectFactory.createTotalChargeType();
                    totalChargeType.setValue(totalCharge);
                    ChargesType chargesType=objectFactory.createChargesType();
                    chargesType.setTotalCharge(totalChargeType);
                    divcon.setCharges(chargesType);
                    divcon.getCmlOrTargetOrLigand().add(cml);
                    JAXBContext jc = JAXBContext.newInstance("com.quantumbioinc.xml.divcon");
                    Marshaller m = jc.createMarshaller();
                    m.setProperty(Marshaller.JAXB_SCHEMA_LOCATION, "http://quantumbioinc.com/schema/divcon /home/roger/NetBeansProjects/OOBackbone/schemas/divcon.xsd");
                    m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
                    DivconNamespacePrefixMapper dnpm=new DivconNamespacePrefixMapper();
                    m.setProperty("com.sun.xml.bind.namespacePrefixMapper", dnpm);
            //Marshal object into file.
                    ByteArrayOutputStream baos=new java.io.ByteArrayOutputStream();
                    JAXBElement<com.quantumbioinc.xml.divcon.DivconType> jaxbOutElement=objectFactory.createDivcon(divcon);
                    m.marshal(jaxbOutElement, baos);
//                    TopologySymmetry topologySymmetry=new TopologySymmetry();
                    ArrayList classNumbers=new ArrayList();//topologySymmetry.perceive(baos.toString());
                    int[] cn=new int[classNumbers.size()];
                    for(int index=0;index<classNumbers.size();index++)
                    {
                        cn[index]=(int)classNumbers.get(index);
                    }
//                    System.out.println("cn: "+cn.length);
        SVLVar[] model=new SVLVar[4];
        model[0]=svlMol.peek(1);
        model[1]=svlMol.peek(2);
        model[2]=svlMol.peek(3);
        SVLVar[] atomVectors=new SVLVar[svlMol.peek(4).length()];
        for(int index=0;index<svlMol.peek(4).length();index++)
        {
            atomVectors[index]=svlMol.peek(4).peek(index+1);
        }
        model[3]=new SVLVar(atomVectors);
        SVLVar[] data=new SVLVar[2];
        data[0]=new SVLVar(model);
        SVLVar[] cnVar=new SVLVar[1];
        cnVar[0]=new SVLVar(cn);
        data[1]=new SVLVar(new String[]{"qbData.classNumber"},cnVar);
        return new SVLVar(new SVLVar(data), new SVLVar("",true));
    }

    DivconType loadTarget(H5File h5File, String target) throws Exception
    {
                HObject obj=h5File.get("DivCon/"+target+"/Documents/"+target+".xml");
//        java.lang.System.out.println("looking for: "+"DivCon/"+target+"/Documents/"+target+".xml");
        SVLVar[] model=new SVLVar[4];
        if(obj!=null)
        {
            H5ScalarDS doc=(H5ScalarDS)obj;
        JAXBContext jc = JAXBContext.newInstance("com.quantumbioinc.xml.divcon");
        Unmarshaller um=jc.createUnmarshaller();
        JAXBElement<com.quantumbioinc.xml.divcon.DivconType> jaxbOutElement=um.unmarshal(new StreamSource(new ByteArrayInputStream (((String[])doc.read())[0].getBytes())), com.quantumbioinc.xml.divcon.DivconType.class);
        return jaxbOutElement.getValue();
        }
        return null;
    }
    
    void getTargetCoordinates(DivconType divcon, ArrayList<Double> xsList,  ArrayList<Double> ysList,  ArrayList<Double> zsList)
    {
        Cml cml=(Cml)divcon.getCmlOrTargetOrLigand().get(0);
        for(int moleculeCount=0;moleculeCount<cml.getAnyCmlOrAnyOrAny().size();moleculeCount++)
        {
            if(((JAXBElement)cml.getAnyCmlOrAnyOrAny().get(moleculeCount)).getValue() instanceof com.quantumbioinc.xml.divcon.Symmetry) continue;
            JAXBElement<com.quantumbioinc.xml.divcon.Molecule> jaxbMoleculeElement=(JAXBElement<com.quantumbioinc.xml.divcon.Molecule>)cml.getAnyCmlOrAnyOrAny().get(moleculeCount);
            Molecule molecule=jaxbMoleculeElement.getValue();
            for(int submoleculeCount=0;submoleculeCount<molecule.getAnyCmlOrAnyOrAny().size();submoleculeCount++)
            {
                if(((JAXBElement)molecule.getAnyCmlOrAnyOrAny().get(submoleculeCount)).getValue() instanceof com.quantumbioinc.xml.divcon.BondArray)continue;
                JAXBElement<com.quantumbioinc.xml.divcon.Molecule> jaxbSubmoleculeElement=(JAXBElement<com.quantumbioinc.xml.divcon.Molecule>)molecule.getAnyCmlOrAnyOrAny().get(submoleculeCount);
                Molecule submolecule=jaxbSubmoleculeElement.getValue();
                JAXBElement<com.quantumbioinc.xml.divcon.AtomArray> jaxbAtomArrayElement=(JAXBElement<com.quantumbioinc.xml.divcon.AtomArray>)submolecule.getAnyCmlOrAnyOrAny().get(0);
                AtomArray atomArray=jaxbAtomArrayElement.getValue();
                for(int atomCount=0;atomCount<atomArray.getAnyCmlOrAnyOrAny().size();atomCount++)
                {
                    JAXBElement<com.quantumbioinc.xml.divcon.Atom> jaxbAtomElement=(JAXBElement<com.quantumbioinc.xml.divcon.Atom>)atomArray.getAnyCmlOrAnyOrAny().get(atomCount);
                    Atom atom=jaxbAtomElement.getValue();
                    JAXBElement<com.quantumbioinc.xml.divcon.AtomType> jaxbAtomTypeElement=(JAXBElement<com.quantumbioinc.xml.divcon.AtomType>)atom.getAnyCmlOrAnyOrAny().get(0);
                    xsList.add(new Double(atom.getX3().doubleValue()));
                    ysList.add(new Double(atom.getY3().doubleValue()));
                    zsList.add(new Double(atom.getZ3().doubleValue()));
                }
            }
        }
        
    }

}
