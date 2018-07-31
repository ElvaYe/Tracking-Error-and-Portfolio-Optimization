/*
Consider n risky assets, use Markoviz mean - variance to calculate global minimum portfolio weight for the benchmark.
Instead of minimizing portfolio variance, find optimial weight for the portfolio return minus the weighted benchmark return
to minimize tracking error.
*/

#pragma warning(disable: 4819)
#include <ql/quantlib.hpp>
#include<iostream>
#include <fstream>
#include <cmath>
#include "CSVParser.hpp"


using namespace QuantLib;
using namespace std;

Real MarkovizWeight(int n, Matrix mu, Matrix u, Matrix sigma, Real mu_p, Matrix &weight_ini, Real &opt_mu, Real &opt_sigmasq) {
	Real sigma_p;

	Real A = (transpose(mu)*inverse(sigma)*mu)[0][0];
	Real B = (transpose(mu)*inverse(sigma)*u)[0][0];
	Real C = (transpose(u)*inverse(sigma)*u)[0][0];

	Real D = A - B * B / C;

	opt_mu = A / C;
	opt_sigmasq = 1 / C;
	weight_ini = inverse(sigma)*u / C * (A - B * mu_p) / D + inverse(sigma)*mu / B * (mu_p*B - B * B / C) / D;
	sigma_p = (mu_p - B / C)*(mu_p - B / C) / D + 1 / C;

	return sigma_p;
}

void someportfolio(int n, Matrix mu, Matrix sigma, Matrix someweight, Real &somemu, Real &somesigma) {
	somemu = (transpose(someweight)*mu)[0][0];
	somesigma = (transpose(someweight)*sigma*someweight)[0][0];
}


int main() {
	/*read return data*/
	Parser portfolio = Parser("12-stock portfolio.csv");
	int startRow = 0;
	int endRow = 847;
	int row = endRow-startRow+1;//how many rows of return data do we have
	int n = 12;// how many stocks do we have in the portfolio
	
	SequenceStatistics ss;
	SequenceStatistics ssd;
	vector<Real> daily_price_index;
	vector<Real> price_diff_index;

	for (int i = startRow; i<=endRow; i++)
	{
		daily_price_index.clear();
		price_diff_index.clear();
		for (int j = 0; j < n; j++)
		{
			daily_price_index.push_back(stod(portfolio[i][j + 2]));
			price_diff_index.push_back(stod(portfolio[i][j + 2]) - stod(portfolio[i][n + 2]));
		}
		ss.add(daily_price_index);
		ssd.add(price_diff_index);
	}
	//write to txt file
	ofstream file;
	//file.open("output TE.csv", ios::out | ios::trunc);
	file.open("output update weights.csv", ios::out | ios::trunc);

	Matrix sigma = ss.covariance();   // covariance matrix for portfoilo return 
	Matrix sigmad = ssd.covariance(); // covariance matrix for portfoilo excess return


	/* Mean of portfolio return and excess return for each stock*/
	Matrix mu(n, 1);
	Matrix mud(n, 1);
	
	for (int i = 0; i < n; i++)
	{
		mu[i][0] = stod(portfolio[startRow][i + 2]);
		mud[i][0] = stod(portfolio[startRow][i + 2]) - stod(portfolio[startRow][n + 2]);
		for (int j = startRow+1; j<= endRow; j++)
		{
			mu[i][0] += stod(portfolio[j][i+2]);
			mud[i][0] += stod(portfolio[j][i + 2]) - stod(portfolio[j][n + 2]);
		}
		mu[i][0] /= row;
		mud[i][0] /= row;
	}
	//cout << mu[n-1][0] << endl << endl;

	/* Mean of SPY daily return */
	Real mu_SPY=0;

	for (int j = startRow; j <= endRow; j++)
	{
		mu_SPY += stod(portfolio[j][n + 2]);
	}
	mu_SPY /= row;

	/* Tracking Error Minimization. Calculate optimal portfolio using markoviz mean-variance framework*/
	Real mu_p = 0.0013; //give some expected return to get minimum variance weight
	Matrix u(n, 1);
	for (int i = 0; i < n; i++) { u[i][0] = 1; }
	Matrix weight_ini(n, 1);
	for (int i = 0; i < n; i++) { weight_ini[i][0] = 1 / n; } // initialize equal weight for each stock
	
	Matrix weight(n, 1);
	Real sigma_port;
	Real sigma_calopt;
	Real opt_sigmasq = 0.001;
	Real opt_mu = 0;

	sigma_port = MarkovizWeight(n, mud, u, sigmad, mu_p, weight_ini, opt_mu, opt_sigmasq);
	cout << "The weight for portfolio is " << endl;
	cout << weight_ini << endl;
	file << "The weight for portfolio is " << endl;
	for (int i = 0; i < n; i++)
	{
		file << weight_ini[i][0] << endl;
	}
	
	cout << "Average Daily Return is: " << "," << (transpose(weight_ini)*mu)[0][0] << endl << "Average Monthly Return is: " << "," << pow((1 + (transpose(weight_ini)*mu)[0][0]), 22) - 1 << endl
		<< "Average SPY Daily Return is: " << "," << mu_SPY << endl << "Average SPY Monthly Return is: " << "," << pow((1 + mu_SPY), 22) - 1 << endl
		<< "Average Daily Volatility is: " << "," << sqrt((transpose(weight_ini)*sigma*weight_ini)[0][0]) << endl << "Average Monthly Volatility is: " << "," << sqrt((transpose(weight_ini)*sigma*weight_ini)[0][0] * 22) << endl
		<< "Tracking Error of portfolio is: " << "," << sqrt((transpose(weight_ini)*sigmad*weight_ini)[0][0]) << endl << endl;
	
	file << "Average Daily Return is: " <<"," << (transpose(weight_ini)*mu)[0][0] << endl << "Average Monthly Return is: " << "," << pow((1+(transpose(weight_ini)*mu)[0][0]),22)-1 << endl
		<< "Average SPY Daily Return is: " << "," << mu_SPY << endl << "Average SPY Monthly Return is: " << "," << pow((1+mu_SPY), 22) - 1 << endl
		<< "Average Daily Volatility is: " << "," << sqrt((transpose(weight_ini)*sigma*weight_ini)[0][0]) << endl << "Average Monthly Volatility is: " << "," << sqrt((transpose(weight_ini)*sigma*weight_ini)[0][0]*22) << endl
		<< "Tracking Error of portfolio is: " << "," << sqrt((transpose(weight_ini)*sigmad*weight_ini)[0][0]) << endl << endl;

	/* General MPT method */
	Matrix weight_ini_mpt(n, 1);
	for (int i = 0; i < n; i++) { weight_ini_mpt[i][0] = 1 / n; }
	Real mu_mpt = mu_p + mu_SPY;
	Real sigma_port_mpt;
	
	sigma_port_mpt = MarkovizWeight(n, mu, u, sigma, mu_mpt, weight_ini_mpt, opt_mu, opt_sigmasq);
	
	cout << "The mpt weight for portfolio is " << endl;
	cout << weight_ini_mpt << endl;
	
	file << "The mpt weight for portfolio is " << endl;
	for (int i = 0; i < n; i++)
	{
		file << weight_ini_mpt[i][0] << endl;
	}
	file << endl;
	
	cout << "Average Daily Return is: " << "," << (transpose(weight_ini_mpt)*mu)[0][0] << endl << "Average Monthly Return is: " << "," << pow((1 + (transpose(weight_ini_mpt)*mu)[0][0]), 22) - 1 << endl
		<< "Average SPY Daily Return is: " << "," << mu_SPY << endl << "Average SPY Monthly Return is: " << "," << pow((1 + mu_SPY), 22) - 1 << endl
		<< "Average Daily Volatility is: " << "," << sqrt((transpose(weight_ini_mpt)*sigma*weight_ini_mpt)[0][0]) << endl << "Average Monthly Volatility is: " << "," << sqrt((transpose(weight_ini_mpt)*sigma*weight_ini_mpt)[0][0] * 22) << endl << endl;
	file << "Average Daily Return is: " << "," << (transpose(weight_ini_mpt)*mu)[0][0] << endl << "Average Monthly Return is: " << "," << pow((1 + (transpose(weight_ini_mpt)*mu)[0][0]), 22) - 1 << endl
		<< "Average SPY Daily Return is: " << "," << mu_SPY << endl << "Average SPY Monthly Return is: " << "," << pow((1 + mu_SPY), 22) - 1 << endl
		<< "Average Daily Volatility is: " << "," << sqrt((transpose(weight_ini_mpt)*sigma*weight_ini_mpt)[0][0]) << endl << "Average Monthly Volatility is: " << "," << sqrt((transpose(weight_ini_mpt)*sigma*weight_ini_mpt)[0][0] * 22) << endl << endl;
		


	/* Generate the efficient frontier */
	int ct;
	ct = 60;

	Matrix exp_returnp(ct, 1);
	exp_returnp[0][0] = -0.001;
	for (int i = 0; i < ct-1; i++) {
		exp_returnp[i + 1][0] = exp_returnp[i][0] + 0.00005;
	}

	Matrix sigmaTE(ct, 1);
	for (int j = 0; j < ct; j++) {
		sigmaTE[j][0] = MarkovizWeight(n, mud, u, sigmad, exp_returnp[j][0], weight_ini, opt_mu, opt_sigmasq);
	}
	Matrix sigmap(ct, 1);
	for (int j = 0; j < ct; j++) {
		sigmap[j][0] = MarkovizWeight(n, mu, u, sigma, exp_returnp[j][0], weight_ini, opt_mu, opt_sigmasq);
	}

	cout << "Expected Return: " << endl;
	cout << exp_returnp << endl;

	file << "Efficient Frontior" << endl <<"Expected Excess Return: " << endl;
	for (int i = 0; i < ct; i++)
	{
		file << exp_returnp[i][0] << endl;
	}
	cout << "Tracking Error: " << endl;
	cout << sigmaTE << endl;

	file << "Tracking Error: " << endl;
	for (int i = 0; i < ct; i++)
	{
		file << sqrt(sigmaTE[i][0]) << endl;
	}

	cout << "Portfolio Volatility: " << endl;
	cout << sigmap << endl;

	file << "Portfolio Volatility: " << endl;
	for (int i = 0; i < ct; i++)
	{
		file << sqrt(sigmap[i][0]) << endl;
	}
	

	file.close();

	system("pause");
	return 0;

}