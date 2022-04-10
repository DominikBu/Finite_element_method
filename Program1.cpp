#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

struct Node
{
	double x;
	double y;
	int BC;
	double t;
};

struct Element
{
	int id[4];
	double H[4][4];
	double Hbc[4][4];
	double P[4];
	double C[4][4];

	Element()
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				H[i][j] = 0.0;
			}
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				Hbc[i][j] = 0.0;
			}
		}

		for (int i = 0; i < 4; i++)
		{
			P[i] = 0.0;
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				C[i][j] = 0.0;
			}
		}
	}

	void zeruj()
	{
		
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				H[i][j] = 0.0;
			}
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				Hbc[i][j] = 0.0;
			}
		}

		for (int i = 0; i < 4; i++)
		{
			P[i] = 0.0;
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				C[i][j] = 0.0;
			}
		}//*/
	}
};

struct Grid
{
	double H;
	double B;
	int nH;
	int nB;
	int nN;				
	int nE;				
	Node * nodes;
	Element * elements;

	Grid(double h, double b, int nh, int nb, double t0)
	{
		H = h;
		B = b;
		nH = nh;
		nB = nb;
		nN = nH * nB;
		nE = (nH - 1)*(nB - 1);

		nodes = new Node[nN];
		elements = new Element[nE];

		double dx = B / (nB - 1);
		double dy = H / (nH - 1);

		//Uzupełnienie węzłów
		for (int i = 0; i < nN; i++)
		{
			nodes[i].t = t0;
			
			nodes[i].x = (i / nH)*dx;
			nodes[i].y = (i % nH)*dy;
			if (nodes[i].x == 0.0 || nodes[i].x == B || nodes[i].y == 0.0 || nodes[i].y == H)
				nodes[i].BC = 1;
			else//*/
				nodes[i].BC = 0;
		}

		//Uzupełnienie elementów
		int k = 0;
		for (int i = 0; i < nE; i++)
		{
			
			k = (i % (nH - 1)) + (i / (nH - 1)) * nH;
			elements[i].id[0] = k;
			elements[i].id[1] = k + nH;
			elements[i].id[2] = elements[i].id[1] + 1;
			elements[i].id[3] = elements[i].id[0] + 1;
			//*/
		}
	}

	void wyswietlN()
	{
		for (int i = 0; i < nN; i++)
		{
			cout << "N"<<i<<"(" << nodes[i].x << ";" << nodes[i].y << ") " << "BC=" << nodes[i].BC << "\n";
		}
	}

	void wyswietlE()
	{
		for (int i = 0; i < nE; i++)
		{
			cout << "E" << i << "(";
			for (int j = 0; j < 3; j++)
			{
				cout << elements[i].id[j] << ",";
			}
			cout << elements[i].id[3];
			cout << ")\n";
		}
	}
};

struct Gauss
{
	double* w;
	double* A;
	int count;
	Gauss(int n)
	{
		count = n;
		switch (n)
		{
		case 2:
			w = new double[2];
			A = new double[2];
			w[0] = -1 / sqrt(3);
			w[1] = -w[0];
			A[0] = 1;
			A[1] = 1;
			break;
		case 3:
			w = new double[3];
			A = new double[3];
			w[0] = -sqrt(3.0 / 5.0);
			w[1] = 0;
			w[2] = -w[0];
			A[0] = 5.0/9.0;
			A[1] = 8.0 / 9.0;
			A[2] = A[0];
			break;
		}
	}

	
};

double f( double x)
{
	return 5 * x * x + 3 * x + 6;
}

double calkowanie1D(Gauss x, int pkt)
{
	pkt++;
	double suma = 0;
	for (int i = 0; i < pkt; i++)
	{

		suma += x.A[i] * f(x.w[i]);

	}
	return suma;
}

double f2D(double x, double y)
{
	return 5 * x * x * y * y + 3 * x * y + 6;
}

double calkowanie2D(Gauss x, int pkt)
{
	pkt++;
	double suma = 0;
	for (int i = 0; i < pkt; i++)
	{
		for (int j = 0; j < pkt; j++)
		{
			suma += x.A[i] * x.A[j] * f2D(x.w[i], x.w[j] );
		}
	}
	return suma;
}

struct Sciana
{
	double** tab;
	double* waga;
	double* ksi;
	double* eta;

};

struct Element4_2D
{
	int count;
	double ** NKsi;
	double ** NEta;
	double ** N;
	Sciana sciany[4];

	Element4_2D(Gauss x)
	{
		count = x.count;

		NKsi = new double* [count * count];
		for (int i = 0; i < (count * count); i++) NKsi[i] = new double[count * 2];

		NEta = new double* [count * count];
		for (int i = 0; i < (count * count); i++) NEta[i] = new double[count * 2];

		double** tabTemp;
		tabTemp = new double* [count * count];
		for (int i = 0; i < (count * count); i++)
		{
			tabTemp[i] = new double[3];
		}


		int k = 0;
		int z = 1;
		for (int i = 0; i < count; i++)
		{
			for (int j = 0; j < count; j++)
			{
				NKsi[i * count + j][0] = -(1.0 / 4.0) * (1.0 - x.w[i]);
				NKsi[i * count + j][1] = (1.0 / 4.0) * (1.0 - x.w[i]);
				NKsi[i * count + j][2] = (1.0 / 4.0) * (1.0 + x.w[i]);
				NKsi[i * count + j][3] = -(1.0 / 4.0) * (1.0 + x.w[i]);

				NEta[i * count + j][0] = -(1.0 / 4.0) * (1.0 - x.w[k]);
				NEta[i * count + j][1] = -(1.0 / 4.0) * (1.0 + x.w[k]);
				NEta[i * count + j][2] = (1.0 / 4.0) * (1.0 + x.w[k]);
				NEta[i * count + j][3] = (1.0 / 4.0) * (1.0 - x.w[k]);


				tabTemp[i * count + j][0] = x.w[i];
				tabTemp[i * count + j][1] = x.w[k];
				tabTemp[i * count + j][2] = x.w[k];

				k += z;
			}
			k -= z;
			z *= -1;
		}


		//Warunki brzegowe sciany
		for (int i = 0; i < 4; i++)
		{
			sciany[i].ksi = new double[count];
			sciany[i].eta = new double[count];
			sciany[i].waga = new double[count];
			sciany[i].tab = new double* [count];
			for (int j = 0; j < count; j++)
				sciany[i].tab[j] = new double[4];

			for (int j = 0; j < count; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					sciany[i].tab[j][k] = 0.0;
				}
			}
		}

		for (int i = 0; i < count; i++)
		{
	
			sciany[0].ksi[i] = x.w[i];
			sciany[0].eta[i] = -1;
			sciany[0].waga[i] = x.A[i];
		}
		for (int i = 0; i < count; i++)
		{

			sciany[1].ksi[i] = 1;
			sciany[1].eta[i] = x.w[i];
			sciany[1].waga[i] = x.A[i];
		}
		for (int i = 0; i < count; i++)
		{

			sciany[2].ksi[i] = x.w[i];
			sciany[2].eta[i] = 1;
			sciany[2].waga[i] = x.A[i];
		}
		for (int i = 0; i < count; i++)
		{

			sciany[3].ksi[i] = -1;
			sciany[3].eta[i] = x.w[i];
			sciany[3].waga[i] = x.A[i];
		}

		for (int i = 0; i < 4; i++)
		{
			//cout << "Sciana " << i << endl;
			for (int j = 0; j < count; j++)
			{
				sciany[i].tab[j][0] = 1.0 / 4.0 * (1.0 - sciany[i].ksi[j])*(1.0 - sciany[i].eta[j]);
				sciany[i].tab[j][1] = 1.0 / 4.0 * (1.0 + sciany[i].ksi[j])*(1.0 - sciany[i].eta[j]);
				sciany[i].tab[j][2] = 1.0 / 4.0 * (1.0 + sciany[i].ksi[j])*(1.0 + sciany[i].eta[j]);
				sciany[i].tab[j][3] = 1.0 / 4.0 * (1.0 - sciany[i].ksi[j])*(1.0 + sciany[i].eta[j]);

				/*
				cout << j+1 << "\t";
				cout << sciany[i].ksi[j] << "\t";
				cout << sciany[i].eta[j] << "\t";
				cout << sciany[i].tab[j][0] << "\t";
				cout << sciany[i].tab[j][1] << "\t";
				cout << sciany[i].tab[j][2] << "\t";
				cout << sciany[i].tab[j][3] << endl;
				//*/
			}
			
		}

		//Wartosci funkcji ksztaltu w punktach calkowania
		N = new double* [count * count];
		for (int i = 0; i < (count * count); i++) N[i] = new double[count * 2];

		k = 0;
		z = 1;
		for (int i = 0; i < count; i++)
		{
			for (int j = 0; j < count; j++)
			{
				N[i * count + j][0] = (1.0 / 4.0) * (1.0 - x.w[k]) * (1.0 - x.w[i]);
				N[i * count + j][1] = (1.0 / 4.0) * (1.0 + x.w[k]) * (1.0 - x.w[i]);
				N[i * count + j][2] = (1.0 / 4.0) * (1.0 + x.w[k]) * (1.0 + x.w[i]);
				N[i * count + j][3] = (1.0 / 4.0) * (1.0 - x.w[k]) * (1.0 + x.w[i]);

				k += z;
			}
			k -= z;
			z *= -1;
		}

		
		
	}

	void pri()
	{
		cout << "dN/dKsi:" << endl;
		for (int i = 0; i < (count * count); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << NKsi[i][j] << "\t";
			}
			cout << endl;
		}

		cout << endl << endl;

		cout << "dN/dEta:" << endl;
		for (int i = 0; i < (count * count); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << NEta[i][j] << "\t";
			}
			cout << endl;
		}
		
	}

	void pri2()
	{
		cout << "N:" << endl;
		for (int i = 0; i < (count * count); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << N[i][j] << "\t";
			}
			cout << endl;
		}

	}
};

struct Jakobian
{
	double * J;
	double * J_inv;
	double detJ;

	Jakobian()
	{
		J = new double[4];
		J_inv = new double[4];

		for (int i = 0; i < 4; i++)
		{
			J[i] = 0;
			J_inv[i] = 0;
		}
	}

	void print()
	{
		cout << "Jakobian wyznacznik: " << detJ << endl;
		cout << J[0] << "\t" << J[1] << endl;
		cout << J[2] << "\t" << J[3] << endl;
		cout << "Jakobian odwrotny" << endl;
		cout << J_inv[0] << "\t" << J_inv[1] << endl;
		cout << J_inv[2] << "\t" << J_inv[3] << endl << endl;
	}
};

void jakobian(int nE, int npc, Jakobian & jakob, Element4_2D Elm, Grid siatka)
{
	for (int i = 0; i < 4; i++)
	{
		jakob.J[0] += Elm.NKsi[npc][i] * siatka.nodes[siatka.elements[nE].id[i]].x;
		jakob.J[1] += Elm.NKsi[npc][i] * siatka.nodes[siatka.elements[nE].id[i]].y;
		jakob.J[2] += Elm.NEta[npc][i] * siatka.nodes[siatka.elements[nE].id[i]].x;
		jakob.J[3] += Elm.NEta[npc][i] * siatka.nodes[siatka.elements[nE].id[i]].y;
	}

	jakob.detJ = (jakob.J[0] * jakob.J[3]) - (jakob.J[1] * jakob.J[2]);

	jakob.J_inv[0] = jakob.J[3];
	jakob.J_inv[1] = -jakob.J[1];
	jakob.J_inv[2] = -jakob.J[2];
	jakob.J_inv[3] = jakob.J[0];


	for (int i = 0; i < 4; i++)
	{
		jakob.J_inv[i] *= (1.0 / jakob.detJ);
	}
}

void NdxNdy(int npc, Jakobian jakob, double** Ndx, double** Ndy, Element4_2D Elm)
{

	for (int i = 0; i < 4; i++)
	{
		Ndx[i][npc] = jakob.J_inv[0] * Elm.NKsi[npc][i] +jakob.J_inv[1] * Elm.NEta[npc][i];
		Ndy[i][npc] = jakob.J_inv[3] * Elm.NEta[npc][i] +jakob.J_inv[2] * Elm.NKsi[npc][i];
	}

	
	
}

double dlugosc(double a, double b) { return (a + b) / 2; }

void Hp(double** Ndx, double** Ndy, Grid siatka, double k, Element& H, int npc, double det, Element& C, double cw, double gestosc, Element4_2D e42d, Gauss gauss)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			H.H[j][i] = k*det*((Ndx[i][npc] * Ndx[j][npc])  + (Ndy[i][npc] * Ndy[j][npc])) * gauss.A[npc / gauss.count] * gauss.A[npc % gauss.count];
			C.H[j][i] = det * cw * gestosc * (e42d.N[npc][i] * e42d.N[npc][j]) * gauss.A[npc / gauss.count] * gauss.A[npc % gauss.count];
		}
		
	}
	/*
	cout << "H pkt calkowania " << npc << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << H.H[j][i] << "  ";

		}

		cout << endl;

	}
	cout << endl;
	//*//*
	cout << "C pkt calkowania " << npc << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << C.H[j][i] << "  ";

		}

		cout << endl;

	}
	cout << endl;
	//*/
}

void H(Element & e, Element* Hpc, int nr, Element * Cpc, Gauss gauss)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < (gauss.count * gauss.count); k++)
			{
				e.H[j][i] += Hpc[k].H[j][i];
				e.C[j][i] += Cpc[k].H[j][i];
			}		
		}
	}

	/*
	cout << "H koncowe elementu " << nr << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << e.H[j][i] << "  ";

		}

		cout << endl;
	}
	cout << endl << endl;
	//*/

	/*
	cout << "C koncowe elementu " << nr << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << e.C[j][i] << "  ";

		}

		cout << endl;
	}
	cout << endl << endl;
	*/
}

void Hbc(double alfa, Element4_2D e, Element& wynik, double detj, int sciana, int nrcalkowania)
{
	//cout << alfa << " " << detj  << " " << e.sciany[sciana].waga[nrcalkowania] << " " << sciana << " " << nrcalkowania << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			
			wynik.H[i][j] = alfa * detj * e.sciany[sciana].waga[nrcalkowania] * e.sciany[sciana].tab[nrcalkowania][j] * e.sciany[sciana].tab[nrcalkowania][i];
			//cout << wynik.H[i][j] << " " << detj << " " << e.sciany[sciana].tab[nrcalkowania][i] << " " << e.sciany[sciana].tab[nrcalkowania][j] << " " << e.sciany[sciana].waga[nrcalkowania] << endl;
		}
	}
}

void wektorP(double k, Element4_2D e, Element& wynik, double detj, int sciana, int nrcalkowania, double temp_otoczenia)
{
	for (int i = 0; i < 4; i++)
	{
		wynik.H[i][0] = k * detj * e.sciany[sciana].waga[nrcalkowania] * e.sciany[sciana].tab[nrcalkowania][i] * temp_otoczenia;
	}

}

bool gauss_uklad(int n, double** AB, double* X)
{

	int i, j, k;
	double m, s;
	double eps = 0.000001;
	// eliminacja współczynników

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// wyliczanie niewiadomych

	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}

void uzupelnij(Grid& siatka)
{

	std::fstream plik;
	//plik.open("Test1_4_4.txt", std::ios::in);
	plik.open("Test4_31_31_trapez.txt", std::ios::in);
	string bufor;
	std::stringstream ss;
	if (plik.good() == true)
	{
		do {
			plik >> bufor;
		} while (bufor.compare("*Node")!=0);

		do {
			int i;
			plik >> bufor;
			if (bufor.compare("*Element,") == 0)
				continue;
			double liczba;

			bufor.erase(bufor.size() - 1, 1);
			i = atoi(bufor.c_str());
			i--;
			//cout << "bufor=" << i << endl;
			plik >> bufor;
			bufor.erase(bufor.size() - 1, 1);
			ss << bufor;
			ss >> liczba;
			ss.clear();
			//liczba += 100;
			siatka.nodes[i].x = liczba;
			//cout << "bufor=" << liczba << endl;
			plik >> bufor;
			ss << bufor;
			ss >> liczba;
			ss.clear();
			//liczba += 100;
			siatka.nodes[i].y = liczba;
			
			//cout << "bufor=" << liczba << endl;
		} while (bufor.compare("*Element,") != 0);
		plik >> bufor;
		do {
			int i;
			plik >> bufor;
			if (bufor.compare("*BC") == 0)
				continue;
			bufor.erase(bufor.size() - 1, 1);
			i = atoi(bufor.c_str());
			i--;
			//cout << "bufor=" << i << endl;
			int n;
			for (int j = 0; j < 3; j++)
			{
				plik >> bufor;
				bufor.erase(bufor.size() - 1, 1);
				n = atoi(bufor.c_str());
				n--;
				siatka.elements[i].id[j] = n;
				//cout << "bufor=" << n << endl;
			}

			plik >> bufor;
			n = atoi(bufor.c_str());
			n--;
			siatka.elements[i].id[3] = n;
			//cout << "bufor=" << n << endl;
			
		}while (bufor.compare("*BC") != 0);

		int i;
		do {
			plik >> bufor;
			if (bufor[bufor.length() - 1] != ',' )
				break;
			bufor.erase(bufor.size() - 1, 1);
			i = atoi(bufor.c_str());
			i--;
			siatka.nodes[i].BC = 1;
		} while (true);

		plik >> bufor;
		i = atoi(bufor.c_str());
		i--;
		siatka.nodes[i].BC = 1;
		plik.close();
	}
	else
	{
		cout << "Nie udalo sie otowrzyc pliku" << endl;
	}


}

int main()
{
	
	double h, b, k, cw, gestosc, t0, alfa, dtau;
	int nH, nB, ilosc_pc, iteracje;
	/*
	cout << "Podaj wysokosc: ";
	cin >> h;
	cout << "Podaj szerokosc: ";
	cin >> b;
	cout << "Podaj ilosc wezlow na wysokosci: ";
	cin >> nH;
	cout << "Podaj ilosc wezlow na szerokosci: ";
	cin >> nB;
	*/
	/*
	double temp_otoczenia[4] = { 50.0,50.0,50.0,50.0 };
	alfa = 300.0; //30;		//
	t0 = 100.0; // 21.0;		//
	h = 0.1;			//0.1
	b = 0.1;			//0.1
	k = 25.0;	// 0.025;		//
	cw = 700.0;	//1000;		//700.0;
	gestosc = 1.2;	//1.2;
	dtau = 50;    //50.0;
	ilosc_pc = 2;
	int iteracje = 10;
	nH = 4;		//4
	nB = 4;	//4
	*/

	double temp_otoczenia[4] = { 1200.0, 1200.0, 1200.0, 1200.0};
	alfa = 300.0; 		
	t0 = 100.0; 	
	h = 0.1;			
	b = 0.1;			
	k = 25.0;	
	cw = 700.0;	
	gestosc = 7800.0;	
	dtau = 50.0;    
	ilosc_pc = 2;
	iteracje = 10;
	nH = 4;
	nB = 4;

	/*
	dtau = 50;
	nH = 4;
	nB = 4;
	*/
	
	

	
	//Dane z prezentacji
	/*
	h = 0.025;
	b = 0.025;
	nH = 2;
	nB = 2;
	*/

	Grid siatka(h,b,nH,nB, t0);


	//uzupelnij(siatka);
	//siatka.wyswietlN();
	//cout << endl;
	//siatka.wyswietlE();
	

	Gauss a(ilosc_pc);
	//Gauss b(3);

	/*
	cout << "wynik calki1: " << calkowanie1D(a, 1) << endl;
	cout << "wynik calki2: " << calkowanie1D(b, 2) << endl;
	cout << "wynik calki3: " << calkowanie2D(a, 1) << endl;
	cout << "wynik calki3: " << calkowanie2D(b, 2) << endl;*/

	Element4_2D t1(a);
	//t1.pri2();
	//cout << endl << endl;

	//Element4_2D t2(b);
	//t2.pri();
	//cout << endl << endl;
	//exit(0);

	double** Hg = new double* [nH * nB];
	for (int i = 0; i < (nH * nB); i++) Hg[i] = new double[nH * nB];
	double** Cg = new double* [nH * nB];
	for (int i = 0; i < (nH * nB); i++) Cg[i] = new double[nH * nB];
	double** HgCg = new double* [nH * nB];
	for (int i = 0; i < (nH * nB); i++) HgCg[i] = new double[nH * nB];
	double* P_agregacja = new double[nH * nB];
	double* PCg = new double[nH * nB];
	double* PgCg = new double[nH * nB];

	for (int czas = 0; czas < iteracje; czas++)
	{
		/*
		cout << "Iteracja " << czas << endl;

		for (int i = 0; i < siatka.nE; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				for (int j2 = 0; j2 < 4; j2++)
				{
					cout << siatka.elements[i].C[j][j2] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}
		cout << "Koniec wypisywania dla iteracji " << czas << endl;
		//*/


		for (int i = 0; i < siatka.nE; i++)
		{
			siatka.elements[i].zeruj();
		}//*/
		

		for (int i = 0; i < (nH * nB); i++)
		{
			for (int j = 0; j < (nH * nB); j++)
			{
				Hg[i][j] = 0.0;
			}
		}

		for (int i = 0; i < (nH * nB); i++)
		{
			for (int j = 0; j < (nH * nB); j++)
			{
				Cg[i][j] = 0.0;
			}
		}

		for (int i = 0; i < (nH * nB); i++)
		{
			for (int j = 0; j < (nH * nB); j++)
			{
				HgCg[i][j] = 0.0;
			}
		}

		for (int j = 0; j < (nH * nB); j++)
		{
			P_agregacja[j] = 0.0;
		}

		for (int j = 0; j < (nH * nB); j++)
		{
			PCg[j] = 0.0;
		}

		for (int j = 0; j < (nH * nB); j++)
		{
			PgCg[j] = 0.0;
		}

		Element tempHbc;
		Element resultHbc[4];
		Element tempP;
		Element resultP[4];

		for (int i = 0; i < siatka.nE; i++)
		{
			
			//Element Hpc[4];
			Element* Cpc = new Element[a.count * a.count];
			Element* Hpc = new Element [a.count*a.count];

			double** dNdx = new double* [4];
			for (int j = 0; j < 4; j++) dNdx[j] = new double[a.count*a.count];

			for (int j = 0; j < 4; j++)
			{
				for (int j2 = 0; j2 < (a.count * a.count); j2++)
				{
					dNdx[j][j2] = 0.0;
				}
			}

			double** dNdy = new double* [4];
			for (int j = 0; j < 4; j++) dNdy[j] = new double[a.count * a.count];

			for (int j = 0; j < 4; j++)
			{
				for (int j2 = 0; j2 < (a.count * a.count); j2++)
				{
					dNdy[j][j2] = 0.0;
				}
			}


			for (int j = 0; j < (t1.count * t1.count); j++)
			{
				Jakobian jakob;
				jakobian(i, j, jakob, t1, siatka);
				//jakob.print();

				NdxNdy(j, jakob, dNdx, dNdy, t1);
				Hp(dNdx, dNdy, siatka, k, Hpc[j], j, jakob.detJ, Cpc[j], cw, gestosc, t1, a);
			}

			H(siatka.elements[i], Hpc, i, Cpc, a);
		
			
			//Gdy warunki brzegowe
			int ktora_sciana[4] = { 0,0,0,0 };
			double detj[4] = { 0.0, 0.0, 0.0, 0.0 };
			if (siatka.nodes[siatka.elements[i].id[0]].BC == 1 && siatka.nodes[siatka.elements[i].id[1]].BC == 1)
			{
				ktora_sciana[0] = 1;
				detj[0] = (sqrt(pow((siatka.nodes[siatka.elements[i].id[1]].x - siatka.nodes[siatka.elements[i].id[0]].x), 2) + pow((siatka.nodes[siatka.elements[i].id[1]].y - siatka.nodes[siatka.elements[i].id[0]].y), 2)))/2;
				//detj[0] = abs(siatka.nodes[siatka.elements[i].id[1]].x - siatka.nodes[siatka.elements[i].id[0]].x) / 2;
			}
			if (siatka.nodes[siatka.elements[i].id[1]].BC == 1 && siatka.nodes[siatka.elements[i].id[2]].BC == 1)
			{
				ktora_sciana[1] = 1;
				detj[1] = (sqrt(pow((siatka.nodes[siatka.elements[i].id[2]].x - siatka.nodes[siatka.elements[i].id[1]].x), 2) + pow((siatka.nodes[siatka.elements[i].id[2]].y - siatka.nodes[siatka.elements[i].id[1]].y), 2))) / 2;
				//detj[1] = abs(siatka.nodes[siatka.elements[i].id[1]].x - siatka.nodes[siatka.elements[i].id[2]].x) / 2;
			}
			if (siatka.nodes[siatka.elements[i].id[2]].BC == 1 && siatka.nodes[siatka.elements[i].id[3]].BC == 1)
			{
				ktora_sciana[2] = 1;
				detj[2] = (sqrt(pow((siatka.nodes[siatka.elements[i].id[2]].x - siatka.nodes[siatka.elements[i].id[3]].x), 2) + pow((siatka.nodes[siatka.elements[i].id[2]].y - siatka.nodes[siatka.elements[i].id[3]].y), 2))) / 2;
				//detj[2] = abs(siatka.nodes[siatka.elements[i].id[2]].x - siatka.nodes[siatka.elements[i].id[3]].x) / 2;
			}
			if (siatka.nodes[siatka.elements[i].id[0]].BC == 1 && siatka.nodes[siatka.elements[i].id[3]].BC == 1)
			{
				ktora_sciana[3] = 1;
				detj[3] = (sqrt(pow((siatka.nodes[siatka.elements[i].id[3]].x - siatka.nodes[siatka.elements[i].id[0]].x), 2) + pow((siatka.nodes[siatka.elements[i].id[3]].y - siatka.nodes[siatka.elements[i].id[0]].y), 2))) / 2;
				//detj[3] = abs(siatka.nodes[siatka.elements[i].id[3]].x - siatka.nodes[siatka.elements[i].id[0]].x) / 2;
			}

			//Zerowanie macierzy sciany
			for (int j = 0; j < 4; j++)
			{
				for (int a1 = 0; a1 < 4; a1++)
				{
					for (int a2 = 0; a2 < 4; a2++)
					{
						resultHbc[j].H[a1][a2] = 0.0;
					}
					resultP[j].H[a1][0] = 0.0;
				}
			}

			for (int j = 0; j < 4; j++)
			{
				if (ktora_sciana[j] == 1)
				{
					for (int l = 0; l < t1.count; l++)
					{

						//Zerowanie macierzy tymczasowej
						for (int a1 = 0; a1 < 4; a1++)
						{
							for (int a2 = 0; a2 < 4; a2++)
							{
								tempHbc.H[a1][a2] = 0.0;
							}
							tempP.H[a1][0] = 0.0;
						}


						Hbc(alfa, t1, tempHbc, detj[j], j, l);
						
						/*
						for (int a1 = 0; a1 < 4; a1++)
						{
							for (int a2 = 0; a2 < 4; a2++)
							{
								cout << tempHbc.H[a1][a2] << " ";
							}
							cout << endl;
						}*/

						wektorP(alfa, t1, tempP, detj[j], j, l, temp_otoczenia[j]);

						for (int a1 = 0; a1 < 4; a1++)
						{
							for (int a2 = 0; a2 < 4; a2++)
							{
								resultHbc[j].H[a1][a2] += tempHbc.H[a1][a2];
							}
							resultP[j].H[a1][0] += tempP.H[a1][0];
						}

						
						
					}
				}
			}
			/*
			//Wypisanie Hbc dla kazdej sciany
			cout << "HBC" << i << endl;
			for (int j = 0; j < 4; j++)
			{
				for (int a1 = 0; a1 < 4; a1++)
				{
					for (int a2 = 0; a2 < 4; a2++)
					{
						cout << resultHbc[j].H[a1][a2] << " ";
					}
					cout << endl;
				}
				cout << endl;
			}
			
			//*/

			/*
			//Wyswietlanie wektora p elementow
			for (int j = 0; j < 4; j++)
			{
				cout << siatka.elements[i].P[j] << "\t";
			}
			cout << endl;
			//*/
			//if(i==3)
			
			//Przypisanie wartosci Hbc elementom
			for (int j = 0; j < 4; j++)
			{
				for (int a1 = 0; a1 < 4; a1++)
				{
					for (int a2 = 0; a2 < 4; a2++)
					{
						siatka.elements[i].Hbc[a1][a2] += resultHbc[j].H[a1][a2];
					}
					siatka.elements[i].P[a1] += resultP[j].H[a1][0];
				}
			}

			



			/*
			//Wypisanie bc dla kazdego elementu
			for (int a1 = 0; a1 < 4; a1++)
			{
				for (int a2 = 0; a2 < 4; a2++)
				{
					cout << siatka.elements[i].Hbc[a1][a2] << " ";
				}
				cout << endl;
			}
			cout << endl << endl;
			*/


			//liczenie Hgcg i Pg koncowe
			for (int j = 0; j < 4; j++)
			{
				for (int l = 0; l < 4; l++)
				{
					Hg[siatka.elements[i].id[j]][siatka.elements[i].id[l]] += siatka.elements[i].H[j][l] + siatka.elements[i].Hbc[j][l];
					Cg[siatka.elements[i].id[j]][siatka.elements[i].id[l]] += siatka.elements[i].C[j][l];
					HgCg[siatka.elements[i].id[j]][siatka.elements[i].id[l]] = (Cg[siatka.elements[i].id[j]][siatka.elements[i].id[l]]) / dtau + Hg[siatka.elements[i].id[j]][siatka.elements[i].id[l]];
				}
				P_agregacja[siatka.elements[i].id[j]] += siatka.elements[i].P[j];
			}
		} //Koniec petli po elementach


		for (int j = 0; j < (nH * nB); j++)
		{
			for (int k = 0; k < (nH * nB); k++)
			{
				PCg[j] += ((Cg[j][k]) / dtau) * siatka.nodes[k].t;
			}
			PgCg[j] = P_agregacja[j] + PCg[j];
		}

		/*
		//Wypisanie Hg koncowe
		for (int j = 0; j < (nH * nB); j++)
		{
			for (int l = 0; l < (nH * nB); l++)
			{
				cout << Hg[j][l] << "  ";
			}
			cout << endl;
		}
		cout << endl << endl;
		//*/

		/*
		//Wypisanie Cg koncowe
		for (int j = 0; j < (nH * nB); j++)
		{
			for (int l = 0; l < (nH * nB); l++)
			{
				cout << Cg[j][l] << "  ";
			}
			cout << endl;
		}
		cout << endl << endl;
		//*/

		/*
		//Wypisanie HgCg koncowe
		for (int j = 0; j < (nH * nB); j++)
		{
			for (int l = 0; l < (nH * nB); l++)
			{
				cout << HgCg[j][l] << "  ";
			}
			cout << endl;
		}
		cout << endl;
		//*/

		/*cout << "to jest p" << endl;
		//Wypisywanie P koncowe
		for (int j = 0; j < (nH * nB); j++)
		{
			cout << P_agregacja[j] << "\t";
		}
		cout << endl;
		//*/
		/*
		//Wypisanie PgCg koncowe
		for (int j = 0; j < (nH * nB); j++)
		{
			cout << PgCg[j] << "\t";
		}
		cout << endl;
		//*/

		double** uklad_rownan = new double* [nH * nB];
		for (int i = 0; i < (nH * nB); i++) uklad_rownan[i] = new double[(nH * nB) + 1];

		for (int i = 0; i < (nH * nB); i++)
		{
			for (int j = 0; j < (nH * nB); j++)
			{
				uklad_rownan[i][j] = HgCg[i][j];
			}
		}
		for (int j = 0; j < (nH * nB); j++)
		{
			uklad_rownan[j][(nH * nB)] = PgCg[j];
		}

		/*
		//Macierz do ukladu rownan gaussa
		for (int i = 0; i < (nH * nB); i++)
		{
			for (int j = 0; j < ((nH * nB)+1); j++)
			{
				cout << uklad_rownan[i][j] << " ";
			}
			cout << endl;
		}
		*/

		double* t = new double[nH * nB];
		if (!gauss_uklad((nH * nB), uklad_rownan, t))
		{
			cout << "DZIELNIK ZERO, ERROR GAUSS\n";
			return 0;
		}


		double min, max;
		min = t[0];
		max = t[0];
		for (int i = 0; i < (nH * nB); i++)
		{
			siatka.nodes[i].t = t[i];
			if (min > t[i])
				min = t[i];
			if (max < t[i])
				max = t[i];
		}

		cout << "iteracja: " << czas << " min: " << min << " max: " << max << endl;

		/*
		//Wartosci temperatur w wezlach
		for (int i = 0; i < (nH * nB); i++)
			cout << siatka.nodes[i].t << endl;

		cout << endl;
		//*/
		
	}
	


	return 0;
}