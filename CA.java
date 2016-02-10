import java.io.*;
import java.util.*;
import static java.lang.Math.*;

public class CA {
	public static String version="8th Feb 2016";
	MersenneTwisterFast random;
	public int mag=1;
	public boolean chemo=false; // By default we are not applying chemo
	public float [][] Oxygen=null;
	public int [][] Receptors=null;
	boolean [][] Cells=null;  // Cells  and the cell type
	boolean [][] Environment=null; // Is the ECM helping with resistance?
	int [][] Vasculature=null; // Initial vasculature as a boolean-like lattice
	Bag cellList=null; // Used in the iterateCells method. Defined here for performance issues
	boolean finishedRun=false;

	int size = 500; // Size of the system lattice
	int timestep=0; // Current Number of timesteps in the simulation
    int births,deaths;

	// Model Parameters
	float initOxygen=0.13f;

	//HERE
	float consumptionBasal=   0.0000f;
	float consumptionDivision=0.0000f; // The Oxygen consumption for cells that proliferate
    float hypoxia=            1.1f; // Threshold
	//HERE

	float pMotility=0.00f; // COMPLETELY ARBITRARY INDEED
	float densityVasculature=0.04f;
	float avgVesselRad=2;

	// Probability of stem cells producing progenitor cells given by oxygen concentration.
	// Probability of progenitor cells of producing differentiated cells: same

	public CA ()
	{
		int time = (int) System.currentTimeMillis();
		random = new MersenneTwisterFast (time);
		//densityVasculature=densityV;
		reset();
		resetVasculature();
	}

	public void reset ()
	{
		Oxygen = new float [size][size];
		Receptors = new int [size][size];
		Environment = new boolean [size][size];
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++) {
				Oxygen[i][j]=initOxygen;
            }

		resetCells();
	}

	void resetVasculature ()
	{
		Vasculature = new int [size][size];

		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++)
				if ((random.nextFloat()<=densityVasculature) && (Vasculature[i][j]==0)) Vasculature[i][j]=1;
	}

	final int distance (int x, int y, int i, int j)
	{
		double dis=Math.sqrt((x-i)*(x-i)+(y-j)*(y-j));

		return (int)Math.round(dis);
	}

	final int[] convertCoordinates(int x, int y)
       {
               // This method converts the coordinates so they take on account
               // the boundaries of the lattice

               if (x < 0) x = size - 1;
               else if (x > size - 1) x = 0;
               if (y < 0) y = size - 1;
               else if (y > size - 1) y = 0;
               int[] result = new int[2];
               result[0] = x; result[1] = y;
               return result;
       }


	public void nextTimeStep ()
	{
		births=0;
		deaths=0;
		//for (int i=0;i<100;i++) iterateOxygen();
		iterateCells();

		//NEW
		int totalCells=0;
		int totalReceptors=0;
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++)
				if (Cells[i][j]) {totalCells++;totalReceptors+=Receptors[i][j];}

        if (totalCells>=1000000) chemo=true;
        timestep++;
        System.out.println ("Cells:\t"+totalCells+" receptors:\t"+((float)totalReceptors/totalCells)+" births:\t"+births+" deaths:\t"+deaths);
	}


	void resetCells ()
	{
		int centre = size/2;
		int radius=50;
		Cells = new boolean [size][size];

		// Let's set the empty space first.
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++) {
				Cells[i][j]=false;
			}

		Cells[centre][centre]=true;
		Receptors[centre][centre]=50;
	}

	public boolean iterateCells()
 	{
		if (cellList==null) cellList = new Bag (size*size);
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++) {
				if (Cells[i][j])  { // All cell types have Cell > 0
					int[] p = new int[2];
					p[0] = i; p[1] = j;
					cellList.add(p);
				}
			}

		while (cellList.size() != 0) {
			// Select the next lattice element at random
			int randomElemIndex=0;
			if (cellList.size()>1) randomElemIndex = random.nextInt(cellList.size()-1);
	   		int[] point = (int[])cellList.get(randomElemIndex);
		   	int rI = point[0];
		   	int rJ = point[1];

			cellList.remove(randomElemIndex); // Remove it from the cell list
		   	//int cell = Cells[rI][rJ];

				if (random.nextFloat()<0.01f) { //10% attrition rate
					Cells[rI][rJ]=false;
					Receptors[rI][rJ]=0;
                    deaths++;
				} else if ((chemo) && (Receptors[rI][rJ]<100*random.nextFloat()) && (random.nextBoolean()) && (!Environment[rI][rJ])) {
					Cells[rI][rJ]=false;
					Receptors[rI][rJ]=0;
				}
			// Do we have space for division or motility ?
			else if (vacantSites(rI,rJ))
				if (Receptors[rI][rJ]<=4*random.nextFloat()*100) {// If tossing the coin we are to proliferate...
						if (Oxygen[rI][rJ]>consumptionDivision) { // AND the oxygen concentration is enough for division...
							Oxygen[rI][rJ]=Oxygen[rI][rJ]-consumptionDivision;
							// Producing daughter cell in an empty neigbhbouring random site
							int[] daughter = findEmptySite (rI,rJ);
							Cells[daughter[0]][daughter[1]]=true;
							Receptors[daughter[0]][daughter[1]]=Receptors[rI][rJ];
							// Mutations? 10%
							if (random.nextFloat()<0.1f) {
								float f = random.nextFloat();
								if (f<0.25f) Receptors[daughter[0]][daughter[1]]--;
								else if (f<0.5f) Receptors[daughter[0]][daughter[1]]++;
								else if (f<0.75f) Receptors[rI][rJ]--;
								else Receptors[rI][rJ]++;
							}
							if (Receptors[daughter[0]][daughter[1]]<=0) Receptors[daughter[0]][daughter[1]]=1;
							else if (Receptors[daughter[0]][daughter[1]]>=100) Receptors[daughter[0]][daughter[1]]=100;
						}
				} else if (pMotility>random.nextFloat()) { // Maybe we can try migration?
						int[] daughter = findEmptySite (rI,rJ);
						Cells[daughter[0]][daughter[1]]=true;
						Cells[rI][rJ]=false;
						Receptors[daughter[0]][daughter[1]]=Receptors[rI][rJ];
						Receptors[rI][rJ]=0;
						System.err.println ("moving "+rI+", "+rJ);
				}
	  	}
 	return true;
 }

final boolean vacantSites (int x, int y)
{
	// Moore neighbourhood
	//int total=0;
	int[] p = new int [2];

	p=convertCoordinates (x+1,y-1);
	if (Cells[p[0]][p[1]]==false) return true;
	p=convertCoordinates (x+1,y);
	if (Cells[p[0]][p[1]]==false) return true;
	p=convertCoordinates (x+1,y+1);
	if (Cells[p[0]][p[1]]==false) return true;
	p=convertCoordinates (x,y-1);
	if (Cells[p[0]][p[1]]==false) return true;
	p=convertCoordinates (x,y+1);
	if (Cells[p[0]][p[1]]==false) return true;
	p=convertCoordinates (x-1,y-1);
	if (Cells[p[0]][p[1]]==false) return true;
	p=convertCoordinates (x-1,y);
	if (Cells[p[0]][p[1]]==false) return true;
	p=convertCoordinates (x-1,y+1);
	if (Cells[p[0]][p[1]]==false) return true;
	return false;
}


int[] findEmptySite (int x, int y)
{
	LinkedList vacantSites = new LinkedList();
	int[] tp1 = new int[2];
	int[] tp2 = new int[2];
	int[] tp3 = new int[2];
	int[] tp4 = new int[2];
	int[] tp5 = new int[2];
	int[] tp6 = new int[2];
	int[] tp7 = new int[2];
	int[] tp8 = new int[2];

	tp1=convertCoordinates (x+1,y-1);
	if (Cells[tp1[0]][tp1[1]]==false) vacantSites.add(tp1);
	tp2=convertCoordinates (x+1,y);
	if (Cells[tp2[0]][tp2[1]]==false) vacantSites.add(tp2);
	tp3=convertCoordinates (x+1,y+1);
	if (Cells[tp3[0]][tp3[1]]==false) vacantSites.add(tp3);
	tp4=convertCoordinates (x,y-1);
	if (Cells[tp4[0]][tp4[1]]==false) vacantSites.add(tp4);
	tp5=convertCoordinates (x,y+1);
	if (Cells[tp5[0]][tp5[1]]==false) vacantSites.add(tp5);
	tp6=convertCoordinates (x-1,y-1);
	if (Cells[tp6[0]][tp6[1]]==false) vacantSites.add(tp6);
	tp7=convertCoordinates (x-1,y);
	if (Cells[tp7[0]][tp7[1]]==false) vacantSites.add(tp7);
	tp8=convertCoordinates (x-1,y+1);
	if (Cells[tp8[0]][tp8[1]]==false) vacantSites.add(tp8);

	// Now let's see where.
	if (vacantSites.size() > 0) { // Now choose a vacant one, otherwise return the original location
		// pick a vacant site and return it
		int vacantElemIndex = random.nextInt(vacantSites.size());
		int[] p = (int[])vacantSites.get(vacantElemIndex);
		return (int[])p;
	} else {
		int[] p = new int[2];
		p[0] = x; p[1] = y; // Just return the original
		System.out.println ("wrong!:"+vacantSites (x,y)+" - "+vacantSites.size());
		return p;
	}

}


/*public void iterateOxygen()
{
	float kDe = 0.001728f;
	float[][] newOxygen = new float[size][size];
	for (int rI = 0; rI < size; rI++)
		for (int rJ = 0; rJ < size; rJ++) {
			// Determine the actual coordinates for top (-1,0), left(0,-1), right(0,1), below(1,0)
			// using periodic boundary conditions
			int[] top = convertCoordinates(rI - 1, rJ);
			int[] left = convertCoordinates(rI, rJ - 1);
			int[] right = convertCoordinates(rI, rJ + 1);
			int[] below = convertCoordinates(rI + 1, rJ);
			// Diffusion
			newOxygen[rI][rJ]
				= Oxygen[rI][rJ] + (kDe *
				(Oxygen[top[0]][top[1]]
				+ Oxygen[left[0]][left[1]]
				+ Oxygen[right[0]][right[1]]
				+ Oxygen[below[0]][below[1]]
				- 4.0f * Oxygen[rI][rJ]));

			// Consumption
			if (Cells[rI][rJ]) newOxygen[rI][rJ]=newOxygen[rI][rJ]-consumptionBasal/10; // Since we iterate 10 times per timestep

			// Production
			if ((Vasculature[rI][rJ]==1) && (Cells[rI][rJ]==false)) {
				newOxygen[rI][rJ]=1.0f;
			}
			// Sanity check
			if (newOxygen[rI][rJ]>1.0f) newOxygen[rI][rJ]=1.0f;
			else if (newOxygen[rI][rJ]<0.0f) newOxygen[rI][rJ]=0.0f;
		}
	Oxygen = newOxygen;
}*/


	public float [][] getOxygen() { return Oxygen; }
	public int [][] getCells () { return Receptors; }
	public int [][] getVasculature() {return Vasculature;}
};
