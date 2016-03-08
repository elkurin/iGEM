/*
 * @file	generation.cpp
 * @brief 	3世代交代ループのシミュレーション
 * 
 * @author 	Eriko Kurimoto <@elkurin>
 * @date 	2016/03/08
 *
 *
 * 詳細説明
 * 世代交代ごとに違う発現をするやつのシミュレーション
 * プロモーター飽和度をn=2, k=0.5のヒル関数で表して、それに比例した発現をしている
 * 3世代すべてにおいて、ANDゲートに必要になるアクチベーターが違うようになっているので、ヤバい
 * すべて同じアクチベーターにするものは、未開発（Aが消滅する前にCが発現してしまい、適切なパラメータがまだ見つかっていない）
 */

#include <iostream>
#include <cmath>
#include <array>
#include <ostream>
#include <fstream>

using namespace std;

namespace {
	ofstream take_log("data_generation.log");
}

//σ因子とアンチσ因子の濃度変数
double A, B, C;
double antiA, antiB, antiC;

//σ因子とアンチ因子の係数
double coefA, coefB, coefC;
double coefAntiA, coefAntiB, coefAntiC;
double coefDecrease;

//分裂時に出てくるアクチベーターの濃度とその係数
double hogeA, hogeB, hogeC;
double increaseHoge;
double decreaseHoge;

//dtの大きさ
const double timeDev = 0.01;

//分裂周期（devideごとに分裂、timeTranの間hogeが発現）
int devide;
int timeTran;

//シミュレーションを走らせる時間
int timeEnd = 1000000;

//何回目のループかのカウンター
int counter;

//アクチベーターヒル関数
double hill1(double k, double b, double x)
{
	return b * x * x / (k * k + x * x);
}

//リプレッサーヒル関数
double hill0(double k, double b, double x)
{
	return b * k * k / (k * k + x * x);
}

//アクチベーターによる飽和度
double f(double i)
{
	double get = hill1(0.5, 1, i);
	// if (get < 0.5) return 0;
	// else return 1;
	return get;
}

//リプレッサーによる飽和度
double g(double i)
{
	double get = hill0(0.5, 1, i);
	// if (get < 0.5) return 0;
	// else return 1;
	return get;
}

//初期値
void init(void)
{
	coefA = 0.15;
	coefB = 0.15;
	coefC = 0.15;

	coefAntiA = 30;
	coefAntiB = 30;
	coefAntiC = 30;

	coefDecrease = 0.1;
	increaseHoge = 1;
	decreaseHoge = 1;

	devide = 100000;
	timeTran = 2500;

	counter = -1;
}

//ある程度以上値小さくなったときに消す丸め込み
void cut()
{
	double min = 10e-8;
	if (A < min) A = 0;
	if (B < min) B = 0;
	if (C < min) C = 0;
	if (antiA < min) antiA = 0;
	if (antiB < min) antiB = 0;
	if (antiC < min) antiC = 0;
	if (hogeA < min) hogeA = 0;
	if (hogeB < min) hogeB = 0;
	if (hogeC < min) hogeC = 0;
}

//微分方程式
void process(void)
{
	double prevA = A, prevB = B, prevC = C;
	double prevAntiA = antiA, prevAntiB = antiB, prevAntiC = antiC;
	double prevHogeA = hogeA, prevHogeB = hogeB, prevHogeC = hogeC;
	hogeA = prevHogeA + (- decreaseHoge * hogeA) * timeDev;
	hogeB = prevHogeB + (- decreaseHoge * hogeB) * timeDev;
	hogeC = prevHogeC + (- decreaseHoge * hogeC) * timeDev;
	A = prevA + (coefA * f(prevA) * g(prevAntiA) - coefDecrease * prevA + coefAntiC * f(prevC) * f(prevHogeC) * g(prevAntiC)) * timeDev;
	B = prevB + (coefB * f(prevB) * g(prevAntiB) - coefDecrease * prevB + coefAntiA * f(prevA) * f(prevHogeA) * g(prevAntiA)) * timeDev;
	C = prevC + (coefC * f(prevC) * g(prevAntiC) - coefDecrease * prevC + coefAntiB * f(prevB) * f(prevHogeB) * g(prevAntiB)) * timeDev;
	antiA = prevAntiA + (coefAntiA * f(prevA) * f(prevHogeA) * g(prevAntiA) - coefDecrease * prevAntiA) * timeDev;
	antiB = prevAntiB + (coefAntiB * f(prevB) * f(prevHogeB) * g(prevAntiB) - coefDecrease * prevAntiB) * timeDev;
	antiC = prevAntiC + (coefAntiC * f(prevC) * f(prevHogeC) * g(prevAntiC) - coefDecrease * prevAntiC) * timeDev;
	
	cut();
}

//メイン関数
int main(void)
{
	init();

	//初期濃度（指定していないものはすべて0）
	A = 0.75;

	//反応ループをtimeEnd * timeDevまで回す
	for (int t = 0; t < timeEnd; t++) {
		if (t % devide == 0) counter++;
		if (t % devide > devide - timeTran) increaseHoge = 1;
		else increaseHoge = 0;
		if (t % devide > devide - timeTran) {
			if (counter % 3 == 0) hogeA += 0.001;
			else if (counter % 3 == 1) hogeB += 0.001;
			else hogeC += 0.001;
		}
		process();
		cout << t * timeDev << " " << A << " " << B << " " << C << " " << antiA << " " << antiB << " " << antiC << " " << hogeA << " " << hogeB << " " << hogeC << endl;
		take_log << A << " " << B << " " << C << " " << antiA << " " << antiB << " " << antiC << " " << hogeA << " " << hogeB << " " << hogeC << endl;
	}

	//gnuplotへのグラフ描画
	FILE *gp;
	gp = popen("gnuplot -persist", "w");
	fprintf(gp, "p \"data_generation.log\" u 0:1 with lines\n"
				"rep \"data_generation.log\" u 0:2 with lines\n"
				"rep \"data_generation.log\" u 0:3 with lines\n"
				"rep \"data_generation.log\" u 0:4 with lines\n"
				"rep \"data_generation.log\" u 0:5 with lines\n"
				"rep \"data_generation.log\" u 0:6 with lines\n"
				"rep \"data_generation.log\" u 0:7 with lines\n");

	pclose(gp);

	return 0;
}

