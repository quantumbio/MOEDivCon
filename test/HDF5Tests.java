/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import svljava.SVLJavaException;
import svljava.SVLVar;

/**
 *
 * @author roger
 */
public class HDF5Tests {
    
         static {
         System.loadLibrary("hdf5");
         }
    public HDF5Tests() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

     @Test
    public void retrieveEigenData() {
         try {
             HDF5Correspondent hdf5Correspondent = new HDF5Correspondent();
             SVLVar cmd = new SVLVar("retrieveEnergyLevels", true);
//             cmd.
             //tokens[1]="s59.pm6.sp";
             SVLVar args = new SVLVar();
             args.=new String[]{"/home/roger/testing-grounds/s59/Moe interface/s59.pm6.sp.xml.h5","s59.pm6.sp"});
             
             SVLVar[] arg = new SVLVar[2];
             arg[0]=cmd;
             arg[1]=new SVLVar(new String[]{"/home/roger/testing-grounds/s59/Moe interface/s59.pm6.sp.xml.h5","s59.pm6.sp"});
             java.lang.System.out.println(arg[1].length()+" "+arg[1].getTokns()[1]);
             int tid = 0;
             SVLVar results = hdf5Correspondent.run(new SVLVar(arg), tid);
         } catch (SVLJavaException sVLJavaException) {
             java.lang.System.out.println(sVLJavaException.getMessage());
         }
     }
}