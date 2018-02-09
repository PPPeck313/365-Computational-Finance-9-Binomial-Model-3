//Preston Peck
//CS 365
//November 20, 2017
//HW9

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;

class Derivative {
    public:
        virtual ~Derivative() {}
        virtual double TerminalPayoff(double S) { return 0; }
        virtual int ValuationTests(double S, double & V) { return 0; }
        virtual void PrintType() {return; }
        // data
        double r;
        double q;
        double sigma;
        double T;

    protected:
        Derivative() { r = 0; q = 0; sigma = 0; T = 0; }
};

class Option : public Derivative {
    public:
        Option() { K = 0; isCall = false; isAmerican = false; }
        virtual ~Option() {}
        virtual double TerminalPayoff(double S);
        virtual int ValuationTests(double S, double &V);
        virtual void PrintType();
        // data
        double K;
        bool isCall;
        bool isAmerican;
};

class BinomialModel {
    public:
        BinomialModel(int n);
        ~BinomialModel();
        int FairValue(int n, Derivative* p_derivative, double S, double t0, double & V);
        int ImpliedVolatility(int n, Derivative* p_derivative, double S, double t0, double target, double& implied_vol, int& num_iter);
    private:
        // methods
        void Clear();
        int Allocate(int n);

        // data
        int n_tree;
        double **stock_nodes;
        double **derivative_nodes;
        int ImpliedVolatilityPrivate(int n, Derivative* p_derivative, double S, double t0, double target, double& implied_vol, int& num_iter);
};

int main() {
    double K, S, sigma, r, T, q = 0.0;
    int n, t0 = 0;

    cout << "S (> 0): $";
    cin >> S;

    while (S <= 0) {
        cout << "S > 0: $";
        cin >> S;
    }

    cout << "K (> 0): $";
    cin >> K;

    while (K <= 0) {
        cout << "K > 0: $";
        cin >> K;
    }

    cout << "r (> 0): ";
    cin >> r;

    while (r <= 0) {
        cout << "r > 0: ";
        cin >> r;
    }

    cout << "q (> 0): ";
    cin >> q;

    while (q <= 0) {
        cout << "q > 0: ";
        cin >> q;
    }

    cout << "sigma (> 0): ";
    cin >> sigma;

    while (sigma <= 0) {
        cout << "sigma > 0: ";
        cin >> sigma;
    }

    cout << "T (> t0): ";
    cin >> T;

    while (T <= t0) {
        cout << "T > 0: ";
        cin >> T;
    } 

    cout << "n (> 0): ";
    cin >> n;

    while (n <= 0) {
        cout << "n > 0: ";
        cin >> n;
    } 

    cout << endl << endl;

    Option Eur_put;
    Eur_put.r = r;
    Eur_put.q = q;
    Eur_put.sigma = sigma;
    Eur_put.T = T;
    Eur_put.K = K;
    Eur_put.isCall = false;
    Eur_put.isAmerican = false;

    Option Eur_call;
    Eur_call.r = r;
    Eur_call.q = q;
    Eur_call.sigma = sigma;
    Eur_call.T = T;
    Eur_call.K = K;
    Eur_call.isCall = true;
    Eur_call.isAmerican = false;

    double FV_Eur_put = 0;
    double FV_Eur_call = 0;
    BinomialModel binom(n);

    double vol = 0.0;
    int iter = 0;
    binom.FairValue(n, &Eur_put, S, t0, FV_Eur_put);
    cout << "European PUT: $" << FV_Eur_put << endl;
    double target = FV_Eur_put + ((S * exp(-q * (T - t0))) - (K * exp(-r * (T - t0))));
    cout << "CALL Target: $" << target << endl << endl;
    binom.ImpliedVolatility(n, &Eur_call, S, t0, target, vol, iter);
    cout << "Sigma = Implied Volatility" << endl;
    cout << sigma << " = " << vol << endl << endl << endl;

    binom.FairValue(n, &Eur_call, S, t0, FV_Eur_call);
    cout << "European CALL: $" << FV_Eur_call << endl;
    target = FV_Eur_call - ((S * exp(-q * (T - t0))) - (K * exp(-r * (T - t0))));
    cout << "PUT Target: $" << target << endl << endl;
    binom.ImpliedVolatility(n, &Eur_put, S, t0, target, vol, iter);
    cout << "Sigma = Implied Volatility" << endl;
    cout << sigma << " = " << vol << endl;
}

double Option::TerminalPayoff(double S) {
    if (isCall) {
        if (S > K) {
            return S - K;
        }
    }
        
    else {
        if (S < K) {
            return K - S;
        }
    }
    return 0;
}

int Option::ValuationTests(double S, double &V) {    
    // early exercise test
    if (isAmerican) {
        if (isCall) {
            V = fmax(V, fmax(S - K, 0));
        }

        else {
            V = fmax(V, fmax(K - S, 0));
        }
    }
    return 0;
}

void Option::PrintType() {
    if (isAmerican) {
        cout << "American ";
    }

    else {
        cout << "European ";
    }

    if (isCall) {
        cout << "CALL";
    }

    else {
        cout << "PUT";
    }
}

int BinomialModel::ImpliedVolatility(int n, Derivative* p_derivative, double S, double t0, double target, double& implied_vol, int& num_iter) {
  int rc = 0;
  const double saved_vol = p_derivative->sigma;
  rc = ImpliedVolatilityPrivate(n, p_derivative, S, t0, target, implied_vol, num_iter);
  p_derivative->sigma = saved_vol;
  return rc;
}

int BinomialModel::ImpliedVolatilityPrivate(int n, Derivative* p_derivative, double S, double t0, double target, double& implied_vol, int& num_iter) {
    cout << "Implied Volatility on ";
    p_derivative->PrintType();
    cout << endl;

    const double tol = pow(1.0, -4);
    const int max_iter = 100;

    implied_vol = 0;
    num_iter = 0;

    double FV = 0.0;
    double sigma_low = 0.01;
    double FV_low = 0.0;
    cout << "Low: $";
    p_derivative->sigma = sigma_low;
    FairValue(n, p_derivative, S, t0, FV_low);
    double diff_FV_low = FV_low - target;
    cout << FV_low << endl;
	
	if (fabs(diff_FV_low) <= tol) {
		implied_vol = p_derivative->sigma;
		cout << endl << "Volatility Rate: " << implied_vol << endl << endl;
		return 0;
	}
	
    double sigma_high = 3.0;
    double FV_high = 0.0;
    cout << "High: $";
    p_derivative->sigma = sigma_high;
	FairValue(n, p_derivative, S, t0, FV_high);
    double diff_FV_high = FV_high - target;
    cout << FV_high << endl;
	
	if (abs(diff_FV_high) <= tol) {
		implied_vol = p_derivative->sigma;
		cout << endl <<  "Volatility Rate: " << implied_vol << endl;
		return 0;
    }
	
	if ((diff_FV_low * diff_FV_high) > 0) {
		implied_vol = 0;
		cout << endl <<  "Volatility Rate: " << implied_vol << endl << endl;
		return 1;
    }
    
    cout << endl;
	
	for (int i = 0; i < max_iter; ++i) {
		p_derivative->sigma = (sigma_low + sigma_high) / 2.0;
        FairValue(n, p_derivative, S, t0, FV);
        double diff_FV = FV - target;
		
		if (abs(diff_FV) <= tol || sigma_high - sigma_low <= tol) {
            implied_vol = p_derivative->sigma;
            num_iter = i;
            cout << endl << "Iterations: " << num_iter << endl;
			cout << "Fair Value: $" << FV << endl;
            cout << "Volatility Rate: " << implied_vol << endl << endl;
            return 0;
		}
		
		else if ((diff_FV_low * diff_FV) > 0) {
            cout << "Fair Value: $" << FV << endl;
            cout << "LOWER: " << p_derivative->sigma << endl;
			sigma_low = p_derivative->sigma;
		}
		
		else {
            cout << "Fair Value: $" << FV << endl;
            cout << "UPPER: " << p_derivative->sigma << endl;
			sigma_high = p_derivative->sigma;
		}
	}
	
    num_iter = max_iter;
    implied_vol = 0;
	cout << endl << endl <<  "Volatility Rate: " << implied_vol << endl << endl;
    return 1;
}

BinomialModel::BinomialModel(int n) {
    n_tree = 0;
    stock_nodes = 0;
    derivative_nodes = 0;
    Allocate(n);
}

BinomialModel::~BinomialModel() {
    Clear(); 
}

void BinomialModel::Clear() {
    if (stock_nodes == NULL && derivative_nodes == NULL) {
        return;
    }

    else {
        for (int i = 0; i <= n_tree; ++i) {
            delete[] stock_nodes[i];
            delete[] derivative_nodes[i];
        }

        delete[] stock_nodes;
        delete[] derivative_nodes;
    }
}

int BinomialModel::Allocate(int n) {
    if (n <= n_tree) return 0;
    // deallocate old tree
    Clear();
    // allocate memory
    n_tree = n;
    stock_nodes = new double* [n + 1];
    derivative_nodes = new double* [n + 1];

    for (int i = 0; i <= n_tree; ++i) {
        stock_nodes[i] = new double[n + 1];
        derivative_nodes[i] = new double[n + 1];

        double* S_tmp = stock_nodes[i];
        double* V_tmp = derivative_nodes[i];

        for (int j = 0; j <= n_tree; ++j) {
            S_tmp[j] = 0;
            V_tmp[j] = 0;
        }
    }
    return 0;
}

int BinomialModel::FairValue(int n, Derivative* p_derivative, double S, double t0, double & V) {     
    int rc = 0;

    V = 0;

    // validation checks
    if (n < 1 || 
        S <= 0 || 
        p_derivative == NULL || 
        p_derivative->T <= t0 ||
        p_derivative->sigma <= 0.0) {
        return 0;
    }
    
    // declaration of local variables (I use S_tmp and V_tmp)
    double* S_tmp = NULL;
    double* V_tmp = NULL;

    double r = p_derivative->r;
    double q = p_derivative->q;
    double T = p_derivative->T;
    double sigma = p_derivative->sigma;

    // calculate parameters
    double dt = (T - t0) / double(n);
    double df = exp(-r * dt);
    double growth = exp((r - q) * dt);
    double u = exp(sigma * sqrt(dt));
    double d = 1.0 / u;
    double p_prob = (growth - d)/(u - d);
    double q_prob = 1.0 - p_prob;

    // more validation checks
    if (p_prob < 0.0 || p_prob > 1.0) {
        return 1;
    }

    // allocate memory if required (call Allocate(n))
    Allocate(n);

    // set up stock prices in tree
    S_tmp = stock_nodes[0];
    S_tmp[0] = S;

    for (int i = 1; i <= n; ++i) {
        double* prev = stock_nodes[i - 1];
        S_tmp = stock_nodes[i];
        S_tmp[0] = prev[0] * d;

        for (int j = 1; j <= n; ++j) {
            S_tmp[j] = S_tmp[j - 1] * u * u;
        }
    }

    // set terminal payoff  (call virtual function in derivative class to calculate payoff)
    int i = n;//n is current scope, n_tree represent maximum tree size
    S_tmp = stock_nodes[i];
    V_tmp = derivative_nodes[i];
    
    for (int j = 0; j <= n; ++j) {
        V_tmp[j] = p_derivative->TerminalPayoff(S_tmp[j]);
    }

    // valuation loop  (call virtual function in derivative class for valuation tests)
    for (int i = n-1; i >= 0; --i) {
        S_tmp = stock_nodes[i];
        V_tmp = derivative_nodes[i];
        double* V_next = derivative_nodes[i + 1];

        for (int j = 0; j <= i; ++j) {
            V_tmp[j] = df * (p_prob * V_next[j + 1] + q_prob * V_next[j]);
            p_derivative->ValuationTests(S_tmp[j], V_tmp[j]);  // VALUATION TESTS
        }
    }
    // option fair value
    V_tmp = derivative_nodes[0];
    V = V_tmp[0];
    return 0; 
}
