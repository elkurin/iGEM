#include <iostream>
#include <cmath>
#include <ostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>

#define rep(i, n) for (int i = 0; i < n; i++)

using namespace std;

random_device rd;

namespace {
	ofstream take_log("data_generation2_cata.log");
	ofstream take_logm("data_generation2m_cata.log");
	ofstream take_coef("data_generation2_coef_cata.log");
}

int numberOfData = 0;
const int runTime = 1;

//諸値
const double dt = 0.0001;
const double cellcycle = 2;
const double cycle = 12;
const int frame = cellcycle * cycle / dt;
const double delay = 0.01;	//転写翻訳にかかる時間
const double delayk = delay / dt;
const double delayr = delay / dt;
const double kanki = 1.2;	//間期の時間の長さ


//濃度
double A[frame], B[frame], C[frame], D[frame], E[frame], F[frame], a[frame], c[frame], e[frame];			//大文字がσ因子、小文字がanti
double mA[frame], mB[frame], mC[frame], mD[frame], mE[frame], mF[frame], ma[frame], mc[frame], me[frame];	//mRNA
double lmA[frame], lmC[frame], lmE[frame];	//DNA複製時に翻訳されるやつlocked
//double amA[frame], amC[frame], amE[frame];	//DNA複製時に翻訳されるやつactivated
double taRNA[frame], mtaRNA[frame];			//lmA, lmC, lmEをactivate
double Pnrd[frame];							//DNA複製時に発現する
double GFP[frame], RFP[frame], mGFP[frame], mRFP[frame];	//タグ

//Group1 -> A, C, E
//Group2 -> B, D, F
//Group3 -> a, c, e

//反応係数
double k1, k2, k3, l, act, kt;
double kGFP, kRFP;

//解離定数
double K1, K2, K3, KP;

//翻訳の速度定数
double r1, r2, r3, rt;
double rGFP, rRFP;

//半減期
double T1, T2, T3;
double Tm1, Tm2, Tm3;
double Tta, Tmta, TPnrd;
double TGFP, TRFP, TmGFP, TmRFP;

double f(double a) 
{
	// double q = 0.5 + (double)(rd() % 1000) * 0.001;
	double q = 1;
	return a * q;
}

//初期化
void init(void)
{
	//初期値
	A[0] = 0.25;
	B[0] = 0;
	C[0] = 0;
	D[0] = 0;
	E[0] = 0;
	F[0] = 0;
	mA[0] = 0;
	mB[0] = 0;
	mC[0] = 0;
	mD[0] = 0;
	mE[0] = 0;
	mF[0] = 0;
	a[0] = 0;
	c[0] = 0;
	e[0] = 0;
	ma[0] = 0;
	mc[0] = 0;
	me[0] = 0;
	lmA[0] = 0;
	lmC[0] = 0;
	lmE[0] = 0;
	// amA[0] = 0;
	// amC[0] = 0;
	// amE[0] = 0;
	taRNA[0] = 0;
	mtaRNA[0] = 0;
	Pnrd[0] = 0;
	GFP[0] = 0;
	RFP[0] = 0;
	mGFP[0] = 0;
	mRFP[0] = 0;

	k1 = f(5);
	k2 = f(5);
	k3 = f(10);
	kt = 10; //固定、未確認
	l = f(9);
	kGFP = 5; //固定
	kRFP = 5; //固定

	K1 = f(1);
	K2 = f(1);
	K3 = f(0.2);

	r1 = f(5);
	r2 = f(5);
	r3 = f(10);
	rt = 10; //固定、未確認
	rGFP = 5; //固定
	rRFP = 5; //固定

	T1 = f(0.13);
	T2 = f(0.13);
	T3 = f(0.13);
	Tm1 = 0.1; //固定
	Tm2 = 0.1; //固定
	Tm3 = 0.1; //固定
	TGFP = 0.2; //固定
	TRFP = 0.2; //固定
	TmGFP = 0.1; //固定
	TmRFP = 0.1; //固定
	Tta = 0.1; //固定、未確認
	Tmta = 0.1; //固定、未確認
	TPnrd = 0.1; //固定、未確認

	take_coef << "No." << numberOfData << " " << k1 << " " << k2 << " " << k3 << " " << l << " " << K1 << " " << K2 << " " << K3 << " " << r1 << " " << r2 << " " << r3 << " " << T1 << " " << T2 << " " << T3 << endl;
}

//自然分解の式
double decrease(double time, double con)
{
	return con - log(2.0) / time * con * dt;
}

void loop(void)
{
	init();

	take_log << A[0] << " " << B[0] << " " << C[0] << " " << D[0] << " " << E[0] << " " << F[0] << " " << a[0] << " " << c[0] << " " << e[0] << " " << GFP[0] << " " << RFP[0] << endl;
	take_logm << mA[0] << " " << mB[0] << " " << mC[0] << " " << mD[0] << " " << mE[0] << " " << mF[0] << " " << ma[0] << " " << mc[0] << " " << me[0] << " " << mGFP[0] << " " << mRFP[0] << " " << lmA[0] << " " << lmC[0] << " " << lmE[0] << endl;

	rep(i, frame - 1) {
		if (i % (frame / 100) == 0) cout << "Simulating : " << i << endl;
		// //前の値を保存しておく
		// double pA = A, pmA = mA;
		// double pB = B, pmB = mB;
		// double pC = C, pmC = mC;
		// double pD = D, pmD = mD;
		// double pE = E, pmE = mE;
		// double pF = F, pmF = mF;
		// double pa = a, pma = ma;
		// double pc = c, pmc = mc;
		// double pe = e, pme = me;
		// double pGFP = GFP, pmGFP = mGFP;
		// double pRFP = RFP, pmRFP = mRFP;

		//自然分解
		A[i + 1] = decrease(T1, A[i]);
		B[i + 1] = decrease(T2, B[i]);
		C[i + 1] = decrease(T1, C[i]);
		D[i + 1] = decrease(T2, D[i]);
		E[i + 1] = decrease(T1, E[i]);
		F[i + 1] = decrease(T2, F[i]);
		a[i + 1] = decrease(T3, a[i]);
		c[i + 1] = decrease(T3, c[i]);
		e[i + 1] = decrease(T3, e[i]);
		mA[i + 1] = decrease(Tm1, mA[i]);
		mB[i + 1] = decrease(Tm2, mB[i]);
		mC[i + 1] = decrease(Tm1, mC[i]);
		mD[i + 1] = decrease(Tm2, mD[i]);
		mE[i + 1] = decrease(Tm1, mE[i]);
		mF[i + 1] = decrease(Tm2, mF[i]);
		ma[i + 1] = decrease(Tm3, ma[i]);
		mc[i + 1] = decrease(Tm3, mc[i]);
		me[i + 1] = decrease(Tm3, me[i]);
		lmA[i + 1] = decrease(Tm1, lmA[i]);
		lmC[i + 1] = decrease(Tm1, lmC[i]);
		lmE[i + 1] = decrease(Tm1, lmE[i]);
		// amA[i + 1] = decrease(Tm1, amA[i]);
		// amC[i + 1] = decrease(Tm1, amC[i]);
		// amE[i + 1] = decrease(Tm1, amE[i]);
		taRNA[i + 1] = decrease(Tta, taRNA[i]);
		mtaRNA[i + 1] = decrease(Tmta, mtaRNA[i]);
		// Pnrd[i + 1] = decrease(TPnrd, Pnrd[i]);
		GFP[i + 1] = decrease(TGFP, GFP[i]);
		RFP[i + 1] = decrease(TRFP, RFP[i]);
		mGFP[i + 1] = decrease(TmGFP, mGFP[i]);
		mRFP[i + 1] = decrease(TmRFP, mRFP[i]);
		
		//転写
		if (i >= delayk) {
			double hilla = a[i - (int)delayk + 1] / (K3 + a[i - (int)delay + 1]);
			double hillc = c[i - (int)delayk + 1] / (K3 + c[i - (int)delay + 1]);
			double hille = e[i - (int)delayk + 1] / (K3 + e[i - (int)delay + 1]);
			double hillA = A[i - (int)delayk + 1] * (1 - hilla) / (K1 + A[i - (int)delayk + 1] * (1 - hilla));
			double hillB = B[i - (int)delayk + 1] / (K2 + B[i - (int)delayk + 1]);
			double hillC = A[i - (int)delayk + 1] * (1 - hillc) / (K1 + C[i - (int)delayk + 1] * (1 - hillc));
			double hillD = D[i - (int)delayk + 1] / (K2 + D[i - (int)delayk + 1]);
			double hillE = E[i - (int)delayk + 1] * (1 - hille) / (K1 + E[i - (int)delayk + 1] * (1 - hille));
			double hillF = F[i - (int)delayk + 1] / (K2 + F[i - (int)delayk + 1]);
			//Pnrdのhillは[Pnrd]が非常に大きいため1とみなせる
			// double hillP = Pnrd[i - (int)delayk + 1] / (KPn + Pnrd[i - (int)delayk + 1]);
			double hillP;
			if (i % (int)(cellcycle / dt) > kanki / dt) {
				//[Pnrd]関連の式が入る
				//瞬時に十分量になって、瞬時に消えるのかよく分からない
				hillP = 1;
			} else hillP = 0;

			mA[i + 1] += k1 * hillA * dt;
			mB[i + 1] += k2 * hillA * dt;
			mC[i + 1] += k1 * hillC * dt;
			mD[i + 1] += k2 * hillC * dt;
			mE[i + 1] += k1 * hillE * dt;
			mF[i + 1] += k2 * hillE * dt;
			ma[i + 1] += k3 * hillC * dt;
			mc[i + 1] += k3 * hillE * dt;
			me[i + 1] += k3 * hillA * dt;
			lmC[i + 1] += l * hillB * dt;
			lmE[i + 1] += l * hillD * dt;
			lmA[i + 1] += l * hillF * dt;
			// amA[i + 1] += act * lmA[i] * taRNA[i] * dt;
			// amC[i + 1] += act * lmC[i] * taRNA[i] * dt;
			// amE[i + 1] += act * lmE[i] * taRNA[i] * dt;
			mtaRNA[i + 1] += kt * hillP * dt;
			mGFP[i + 1] += kGFP * hillB * dt;
			mRFP[i + 1] += kRFP * hillD * dt;
		}

		//翻訳
		if (i >= (int)delayr) {
			A[i + 1] += r1 * mA[i - (int)delayr + 1] * dt;
			B[i + 1] += r2 * mB[i - (int)delayr + 1] * dt;
			C[i + 1] += r1 * mC[i - (int)delayr + 1] * dt;
			D[i + 1] += r2 * mD[i - (int)delayr + 1] * dt;
			E[i + 1] += r1 * mE[i - (int)delayr + 1] * dt;
			F[i + 1] += r2 * mF[i - (int)delayr + 1] * dt;
			a[i + 1] += r3 * ma[i - (int)delayr + 1] * dt;
			c[i + 1] += r3 * mc[i - (int)delayr + 1] * dt;
			e[i + 1] += r3 * me[i - (int)delayr + 1] * dt;
			// A[i + 1] += r1 * amA[i - (int)delayr + 1] * dt;
			// C[i + 1] += r1 * amC[i - (int)delayr + 1] * dt;
			// E[i + 1] += r1 * amE[i - (int)delayr + 1] * dt;
			A[i + 1] += r1 * lmA[i - (int)delayr + 1] * taRNA[i - (int)delayr + 1] * dt;
			C[i + 1] += r1 * lmC[i - (int)delayr + 1] * taRNA[i - (int)delayr + 1] * dt;
			E[i + 1] += r1 * lmE[i - (int)delayr + 1] * taRNA[i - (int)delayr + 1] * dt;
			taRNA[i + 1] += rt * mtaRNA[i - (int)delayr + 1] * dt;
			GFP[i + 1] += rGFP * mGFP[i - (int)delayr + 1] * dt;
			RFP[i + 1] += rRFP * mRFP[i - (int)delayr + 1] * dt;

			if (i % (int)(cellcycle / dt) > kanki / dt) {
				// A[i + 1] += r1 * lmA[i - (int)delayr + 1] * dt;
				// C[i + 1] += r1 * lmC[i - (int)delayr + 1] * dt;
				// E[i + 1] += r1 * lmE[i - (int)delayr + 1] * dt;
			}
		}

		//記録
		take_log << A[i + 1] << " " << B[i + 1] << " " << C[i + 1] << " " << D[i + 1] << " " << E[i + 1] << " " << F[i + 1] << " " << a[i + 1] << " " << c[i + 1] << " " << e[i + 1] << " " << GFP[i + 1] << " " << RFP[i + 1] << endl;
		take_logm << mA[i + 1] << " " << mB[i + 1] << " " << mC[i + 1] << " " << mD[i + 1] << " " << mE[i + 1] << " " << mF[i + 1] << " " << ma[i + 1] << " " << mc[i + 1] << " " << me[i + 1] << " " << mGFP[i + 1] << " " << mRFP[i + 1] << " " << lmA[i + 1] << " " << lmC[i + 1] << " " << lmE[i + 1] << endl;

	}

	//gnuplot出力
	cout << "Outputing to Gnuplot" << endl;
	FILE *gp;
	gp = popen("gnuplot -persist", "w");
	fprintf(gp, 
				// "set terminal postscript color\n"
				// "set output sprintf(\"generetion2_cata_no%d.ps\")\n"
				"set xrange[%d:%d]\n"
				"k = %d\n"
				"p \"data_generation2_cata.log\" u (($0) - k):1 with lines t \"A\"\n"
			    "rep \"data_generation2_cata.log\" u (($0) - k):2 with lines t \"B\"\n"
			    "rep \"data_generation2_cata.log\" u (($0) - k):3 with lines t \"C\"\n"
			    "rep \"data_generation2_cata.log\" u (($0) - k):4 with lines t \"D\"\n"
			    "rep \"data_generation2_cata.log\" u (($0) - k):5 with lines t \"E\"\n"
			    "rep \"data_generation2_cata.log\" u (($0) - k):6 with lines t \"F\"\n", /*numberOfData,*/ 0, frame, (numberOfData - 1) * frame);
	fprintf(gp, 
				"set terminal postscript color\n"
				"set output sprintf(\"generation2_cata_GFP_no%d.ps\")\n"
				"set xrange[%d:%d]\n"
				"k = %d\n"
				"p \"data_generation2_cata.log\" u (($0) - k):10 with lines t \"GFP\"\n"
			    "rep \"data_generation2_cata.log\" u (($0) - k):11  with lines t \"RFP\"\n", numberOfData, 0, frame, (numberOfData - 1) * frame);
	pclose(gp);
}

int main(void)
{
	rep(i, runTime)
	{
		numberOfData++;
		loop();
	}
	return 0;
}
