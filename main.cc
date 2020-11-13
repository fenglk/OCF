/// J. Chem. Phys. 133, 154103 (2010)

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "parameters.h"
#include "FCorrelatorSAt.h"
//#include "FCorrelatorCt.h"
//#include "FCorrelatorCtls.h"


using namespace SAt;
//using namespace Ct;
//using namespace Ctls;


//int main(int argc, char *argv[])
//int main()
int main(int argc, char *argv[])
{
	//string str = "D:\DATA\Gt\100\test\dump\1unwrap.xtc"; /// input xtc file
	//std::string str = "../1unwrap.xtc";
	//string str = "D:\\VS2017\\GtMultiTauCC\\GtMultiTauCC\\munwrap.xtc";
	//char *p = new char[1024];
	//strcpy(p, str.c_str());
	//char *test = "1unwrap.xtc";  /// input xtc file

	XDRFILE *xd;
	int result;
	unsigned int i, j, k;
	
	int natoms;  //total number of atoms
	int Step;
	float Time;
	matrix box;
	//rvec region;
	rvec *x;
	float prec;
	unsigned int counter;
	unsigned int countMol;

	std::cout << argv[0] << "number of argument\t" << argc << std::endl;
	unsigned int nChain11 = atoi(argv[1]);
	unsigned int chainLen11 = atoi(argv[2]);
	unsigned int nChain22 = atoi(argv[3]);
	unsigned int chainLen22 = atoi(argv[4]);

	double nBond1, nBond2, nBond;
	nBond1 = (((double)chainLen11 - 1.0)*((double)nChain11));
	nBond2 = (((double)chainLen22 - 1.0)*((double)nChain22));
	nBond = nBond1 + nBond2;

	int TOTatoms = nChain11 * chainLen11 + nChain22 * chainLen22;
	double phi1, phi2; /// volume fraction for 1 and 2 components
	phi1 = (((double)nChain11)*((double)chainLen11)) / ((double)TOTatoms);
	phi2 =  (((double)nChain22)*((double)chainLen22)) / ((double)TOTatoms);
	unsigned int TOTframes = run / dump;

	/*-------------------------------------*/
	result = read_xtc_natoms(argv[5], &natoms);
	xd = xdrfile_open(argv[5], "r");

	if (exdrOK != result)
	{
		std::cerr << "die_r: read_xtc_natoms\t" << "result:" << result << std::endl;
		exit(1);
	}
	if (NULL == xd)
	{
		std::cerr << "die: xdrfile_open" << std::endl;
		exit(1);
	}
	if (TOTatoms != natoms)
	{
		std::cerr << "Error: incompatile arguments!" << std::endl;
		exit(0);
	}
	//std::cout << "222natoms: " << natoms << std::endl;
	/*-------------------------------------*/


	/*--------------------------------------------------*
	if (exdrOK != result)
		die_r("read_xtc_natoms", result);
	if (NULL == xd)
		die("Opening xdrfile for reading...");
	string readXTCnatoms = "read_xtc_natoms";
	char *chreadXTCnatoms = new char[1024];
	strcpy(chreadXTCnatoms, readXTCnatoms.c_str());
	if (exdrOK != result)
		die_r(chreadXTCnatoms, result);
	/*--------------------------------------------------*/

	counter = 0; //the number of frame reading
	//x = calloc(natoms, sizeof(*x)); //record the coordinations of the atoms in a frame
	x = new rvec[natoms];
	//result = read_xtc(xd, natoms, &Step, &Time, box, x, &prec); //read the number 0 frame from the *.xtc file

	std::ofstream out("readXTClog.txt");
	out << "counter\tStep\tTime\tprec" << std::endl;
	//out << counter << "\t" << Step << "\t" << Time << "\t" << prec << std::endl;

	SAt::Correlator Sxy, Sxz, Syz;
	Sxy.setsize(32, 16, 2, 1); /// only 1 input argument per frame
	Sxy.initialize();
	Sxz.setsize(32, 16, 2, 1);
	Sxz.initialize();
	Syz.setsize(32, 16, 2, 1);
	Syz.initialize();

	SAt::Correlator Sxy1, Sxz1, Syz1;
	Sxy1.setsize(32, 16, 2, 1); /// only 1 input argument per frame
	Sxy1.initialize();
	Sxz1.setsize(32, 16, 2, 1);
	Sxz1.initialize();
	Syz1.setsize(32, 16, 2, 1);
	Syz1.initialize();

	SAt::Correlator Sxy2, Sxz2, Syz2;
	Sxy2.setsize(32, 16, 2, 1); /// only 1 input argument per frame
	Sxy2.initialize();
	Sxz2.setsize(32, 16, 2, 1);
	Sxz2.initialize();
	Syz2.setsize(32, 16, 2, 1);
	Syz2.initialize();


	SAt::Correlator Axy1, Axz1, Ayz1;
	Axy1.setsize(32, 16, 2, nChain11); /// # of chains = input argument per frame
	Axy1.initialize();
	Axz1.setsize(32, 16, 2, nChain11);
	Axz1.initialize();
	Ayz1.setsize(32, 16, 2, nChain11);
	Ayz1.initialize();

	SAt::Correlator Axy2, Axz2, Ayz2;
	Axy2.setsize(32, 16, 2, nChain22); /// # of chains = input argument per frame
	Axy2.initialize();
	Axz2.setsize(32, 16, 2, nChain22);
	Axz2.initialize();
	Ayz2.setsize(32, 16, 2, nChain22);
	Ayz2.initialize();



	/// Orientation Tensor of Chain1 and Chain2: 6 -> xx, yy, zz, xy, xz, yz
	double **OTC1, **OTC2;
	double **OTC12;
	OTC1 = new double*[6];
	OTC2 = new double*[6];
	OTC12 = new double*[6];
	for (i = 0; i < 6; i++) {
		OTC1[i] = new double[nChain11];
		for (j = 0; j < nChain11; j++) {
			OTC1[i][j] = 0.0;
		}
		OTC2[i] = new double[nChain22];
		for (j = 0; j < nChain22; j++) {
			OTC2[i][j] = 0.0;
		}
		OTC12[i] = new double[nChain11 + nChain22];
		for (j = 0; j < (nChain11 + nChain22); j++) {
			OTC12[i][j] = 0.0;
		}
	}


	/// Summation of the Orientation Tensor of Chain1 & Chain2 & (Chain1+Chain2)
	double **SumOTC1, **SumOTC2;
	double **SumOTC;
	SumOTC1 = new double*[6];
	SumOTC2 = new double*[6];
	SumOTC = new double*[6];
	for (i = 0; i < 6; i++) {
		SumOTC1[i] = new double[1];
		SumOTC2[i] = new double[1];
		SumOTC[i] = new double[1];
		SumOTC1[i][0] = 0.0;
		SumOTC2[i][0] = 0.0;
		SumOTC[i][0] = 0.0;
	}





	VecR BondVector;
	unsigned int i1, i2;

	std::ofstream bonglength1("BondLength1.txt");
	std::ofstream bonglength2("BondLength2.txt");
	double length;

	for (i = 0; i < TOTframes + 1; ++i)
	//for (i = 0; i < 1000 + 1; ++i)
	{
		++counter;
		//read a frame from the *.xtc file, x[i][j] is the value of the coordination of the atom i, where j=0,1,2, refers to x,y,z.
		result = read_xtc(xd, natoms, &Step, &Time, box, x, &prec);
		out << counter << "\t" << Step << "\t" << Time << "\t" << prec << std::endl;

		if (i%skipped == 0) {
			std::cout << "current frame\t" << i << std::endl;
			countMol = 0;
			/*-MSD--------------------------------------------------------------------------------------*
			for (j = 0; j < nChain2; ++j) {
				for (k = 0; k < chainLen2; ++k) {
					if ((k >= (chainLen2 - MID) / 2) && (k < (chainLen2 + MID) / 2)) {
						indexAtom = k + j * chainLen2 + nChain1 * chainLen1;
						VSet(molM[countMol], x[indexAtom][0], x[indexAtom][1], x[indexAtom][2]);
						++countMol;
					}
				}
			}
			/*-MSD--------------------------------------------------------------------------------------*/

			/*-EEacf---------------------------------------------------------------------------*
			for (j = 0; j < nChain2; ++j) {
				indexAtom = j * chainLen2 + nChain1 * chainLen1;
				VSet(molM[countMol], x[indexAtom][0] - x[indexAtom + chainLen2 - 1][0], \
												      x[indexAtom][1] - x[indexAtom + chainLen2 - 1][1], \
										             x[indexAtom][2] - x[indexAtom + chainLen2 - 1][2]);
				++countMol;
			}
			/*-EEacf---------------------------------------------------------------------------*/

			/*-Orientation Tensor------------------------------------------------------*/
			for (j = 0; j < 6; j++) {
				for (k = 0; k < nChain11; k++) {
					OTC1[j][k] = 0.0;
				}
			}

			length = 0.0;
			for (j = 0; j < nChain11; j++) {
				for (k = 0; k < chainLen11 - 1; k++) {
					i1 = k + j * chainLen11;
					i2 = i1 + 1;
					VSet(BondVector, x[i1][0] - x[i2][0], x[i1][1] - x[i2][1], x[i1][2] - x[i2][2]); /// Without Normalization by the bond length
					OTC1[0][j] += BondVector.x*BondVector.x;
					OTC1[1][j] += BondVector.y*BondVector.y;
					OTC1[2][j] += BondVector.z*BondVector.z;
					OTC1[3][j] += BondVector.x*BondVector.y;
					OTC1[4][j] += BondVector.x*BondVector.z;
					OTC1[5][j] += BondVector.y*BondVector.z;

					length += VLen(BondVector);		
				}
			}
			bonglength1 << length / nBond1 << std::endl;

			for (j = 0; j < 6; j++) {
				for (k = 0; k < nChain22; k++) {
					OTC2[j][k] = 0.0;
				}
			}

			length = 0.0;
			for (j = 0; j < nChain22; j++) {
				for (k = 0; k < chainLen22 - 1; k++) {
					i1 = k + j * chainLen22 + nChain11 * chainLen11;
					i2 = i1 + 1;
					VSet(BondVector, x[i1][0] - x[i2][0], x[i1][1] - x[i2][1], x[i1][2] - x[i2][2]);
					OTC2[0][j] += BondVector.x*BondVector.x;
					OTC2[1][j] += BondVector.y*BondVector.y;
					OTC2[2][j] += BondVector.z*BondVector.z;
					OTC2[3][j] += BondVector.x*BondVector.y;
					OTC2[4][j] += BondVector.x*BondVector.z;
					OTC2[5][j] += BondVector.y*BondVector.z;

					length += VLen(BondVector);
				}
			}
			bonglength2 << length / nBond2 << std::endl;


			for (k = 0; k < 6; k++) {
				for (j = 0; j < nChain11; j++) {
					OTC12[k][j] = OTC1[k][j];
				}
				for (j = nChain11; j < (nChain11 + nChain22); j++) {
					OTC12[k][j] = OTC2[k][j - nChain11];
				}
			}
		
			Axy1.add(OTC1[3]);
			Axz1.add(OTC1[4]);
			Ayz1.add(OTC1[5]);
			Axy2.add(OTC2[3]);
			Axz2.add(OTC2[4]);
			Ayz2.add(OTC2[5]);

		

			for (j = 0; j < 6; j++) {
				SumOTC1[j][0] = 0.0;
				SumOTC2[j][0] = 0.0;
				SumOTC[j][0] = 0.0;
			}
			
			for (j = 0; j < 6; j++) {
				for (k = 0; k < nChain11; k++) {
					SumOTC1[j][0] += OTC1[j][k];
					SumOTC[j][0] += OTC1[j][k];
				}
				for (k = 0; k < nChain22; k++) {
					SumOTC2[j][0] += OTC2[j][k];
					SumOTC[j][0] += OTC2[j][k];
				}
			}
			/*-Orientation Tensor------------------------------------------------------*/
			
			Sxy.add(SumOTC[3]);
			Sxz.add(SumOTC[4]);
			Syz.add(SumOTC[5]);

			Sxy1.add(SumOTC1[3]);
			Sxz1.add(SumOTC1[4]);
			Syz1.add(SumOTC1[5]);

			Sxy2.add(SumOTC2[3]);
			Sxz2.add(SumOTC2[4]);
			Syz2.add(SumOTC2[5]);
		}
	}

	xdrfile_close(xd);
	out.close();



	Sxy.evaluate();
	Sxz.evaluate();
	Syz.evaluate();
	Sxy1.evaluate();
	Sxz1.evaluate();
	Syz1.evaluate();
	Sxy2.evaluate();
	Sxz2.evaluate();
	Syz2.evaluate();

	Axy1.evaluate();
	Axz1.evaluate();
	Ayz1.evaluate();
	Axy2.evaluate();
	Axz2.evaluate();
	Ayz2.evaluate();

	std::ofstream outall("result_OCFall");
	outall<< "t\t" << "step\t" << "Sxy\tSxz\tSyz\tSxy1\tSxz1\tSyz1\tSxy2\tSxz2\tSyz2\tAxy1\tAxz1\tAyz1\tAxy2\tAxz2\tAyz2\tCxy1\tCxz1\tCyz1\tCxy2\tCxz2\tCyz2\tCxy12\tCxz12\tCyz12\tKappa" << std::endl;
	std::ofstream out1("result_OCFsum");
	out1 << "t\t" << "step\t" << "Stot\tS1\tS2\tA1\tA2\tC1\tC2\tC12\tKappa" << std::endl;
	double st, s1, s2, a1, a2, c1, c2, c12, kappa;

	for (i = 0; i < Sxy.npcorr; ++i)
	{
		st = (Sxy.f[i] + Sxz.f[i] + Syz.f[i]) / nBond / 3.0;
		s1= (Sxy1.f[i] + Sxz1.f[i] + Syz1.f[i]) / nBond1 / 3.0;
		s2 = (Sxy2.f[i] + Sxz2.f[i] + Syz2.f[i]) / nBond2 / 3.0;
		a1 = (Axy1.f[i] + Axz1.f[i] + Ayz1.f[i]) / ((double)chainLen11 - 1.0) / 3.0;
		a2 = (Axy2.f[i] + Axz2.f[i] + Ayz2.f[i]) / ((double)chainLen22 - 1.0) / 3.0;
		c1 = (s1 - a1) / phi1;
		c2 = (s2 - a2) / phi2;
		c12 = (st - phi1 * s1 - phi2 * s2) / (2.0*phi1*phi2);
		kappa = (st - phi1 * a1 - phi2 * a2) / st;
		out1 << Sxy.t[i] * dump*skipped*timestep << "\t" << Sxy.t[i] << "\t"
			<< st  << "\t" << s1 << "\t" << s2 << "\t"
			<< a1 << "\t" << a2 << "\t"
			<< c1 << "\t" << c2 << "\t" << c12 << "\t"
			<< kappa << "\t"
			<< std::endl;


		outall << Sxy.t[i] * dump*skipped*timestep << "\t" << Sxy.t[i] << "\t"
			<< Sxy.f[i] / nBond << "\t"
			<< Sxz.f[i] / nBond << "\t"
			<< Syz.f[i] / nBond << "\t"
			<< Sxy1.f[i] / nBond1 << "\t"
			<< Sxz1.f[i] / nBond1 << "\t"
			<< Syz1.f[i] / nBond1 << "\t"
			<< Sxy2.f[i] / nBond2 << "\t"
			<< Sxz2.f[i] / nBond2 << "\t"
			<< Syz2.f[i] / nBond2 << "\t"
			<< Axy1.f[i] / ((double)chainLen11 - 1.0) << "\t"
			<< Axz1.f[i] / ((double)chainLen11 - 1.0) << "\t"
			<< Ayz1.f[i] / ((double)chainLen11 - 1.0) << "\t"
			<< Axy2.f[i] / ((double)chainLen22 - 1.0) << "\t"
			<< Axz2.f[i] / ((double)chainLen22 - 1.0) << "\t"
			<< Ayz2.f[i] / ((double)chainLen22 - 1.0) << "\t"
			<< (Sxy1.f[i] / nBond1 - Axy1.f[i] / ((double)chainLen11 - 1.0)) / phi1 <<"\t"
			<< (Sxz1.f[i] / nBond1 - Axz1.f[i] / ((double)chainLen11 - 1.0)) / phi1 << "\t"
			<< (Syz1.f[i] / nBond1 - Ayz1.f[i] / ((double)chainLen11 - 1.0)) / phi1 << "\t"
			<< (Sxy2.f[i] / nBond2 - Axy2.f[i] / ((double)chainLen22 - 1.0)) / phi2 << "\t"
			<< (Sxz2.f[i] / nBond2 - Axz2.f[i] / ((double)chainLen22 - 1.0)) / phi2 << "\t"
			<< (Syz2.f[i] / nBond2 - Ayz2.f[i] / ((double)chainLen22 - 1.0)) / phi2 << "\t"
			<< (Sxy.f[i] / nBond - phi1 * Sxy1.f[i] / nBond1 - phi2 * Sxy2.f[i] / nBond2) / (2.0*phi1*phi2) << "\t"
			<< (Sxz.f[i] / nBond - phi1 * Sxz1.f[i] / nBond1 - phi2 * Sxz2.f[i] / nBond2) / (2.0*phi1*phi2) << "\t"
			<< (Syz.f[i] / nBond - phi1 * Syz1.f[i] / nBond1 - phi2 * Syz2.f[i] / nBond2) / (2.0*phi1*phi2) << "\t"
			<< (Sxy.f[i] / nBond - phi1 * Axy1.f[i] / ((double)chainLen11 - 1.0) - phi2 * Axy2.f[i] / ((double)chainLen22 - 1.0)) / (Sxy.f[i] / nBond) << "\t"
			<< (Sxz.f[i] / nBond - phi1 * Axz1.f[i] / ((double)chainLen11 - 1.0) - phi2 * Axz2.f[i] / ((double)chainLen22 - 1.0)) / (Sxz.f[i] / nBond) << "\t"
			<< (Syz.f[i] / nBond - phi1 * Ayz1.f[i] / ((double)chainLen11 - 1.0) - phi2 * Ayz2.f[i] / ((double)chainLen22 - 1.0)) / (Syz.f[i] / nBond) << "\t"
			<< std::endl;


	}

	outall.close();
	out1.close();



	return 0;
}
