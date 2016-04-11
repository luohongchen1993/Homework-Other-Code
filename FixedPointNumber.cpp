#include <iostream>
#include <vector>
#include <string>
using namespace std;

class FixedPointNum {
	// use a static variable to control precision
public:
	FixedPointNum();
	FixedPointNum(const string & numstr);
	FixedPointNum(FixedPointNum & f);
	FixedPointNum(vector<int> n, int s, int ss);

	FixedPointNum& operator=(FixedPointNum &f);
	FixedPointNum& operator+=(FixedPointNum &f);
	FixedPointNum operator+(FixedPointNum &f);
	FixedPointNum& operator-=(FixedPointNum &f);
	FixedPointNum operator-(FixedPointNum &f);
	FixedPointNum& operator*=(FixedPointNum &f);
	FixedPointNum operator*(FixedPointNum &f);
	FixedPointNum& operator/=(FixedPointNum &f);
	FixedPointNum operator/(FixedPointNum &f);

	static void set_precision(int p) { precision = p; };
	static int get_precision() { return precision; };
	void times10();
	void divide10();

	void print();
	int get_scale();
	void change_scale(int s);
	int get_sign();
	vector<int> get_num();
	//friend FixedPointNum add(FixedPointNum f1, FixedPointNum f2); // this function ignores signs
	//FixedPointNum operator+(const FixedPointNum &f); // define this using +=
	// divide, while dec.size()!= precision, *10/divider +mod
//private:
private:
	vector<int> num;
	int scale;
	int sign;
	static int precision;
};
int FixedPointNum::precision = 25;

FixedPointNum::FixedPointNum()
{
	vector<int> new_num;
	new_num.push_back(0);
	num = new_num;
	scale = 0;
	sign = 1;
}

FixedPointNum::FixedPointNum(const string & numstr)
{
	if (numstr[0] == '-')
	{
		sign = -1;
		bool point = false;
		for (int i = 1; i < numstr.size(); i++) {
			if (numstr[i] == '.' && point == false) {
				point = true;
				scale = numstr.size() - i - 1;
				continue;
			}
			else if (isdigit(numstr[i])) {
				int tmp = numstr[i] - '0';
				num.push_back(tmp);
			}
			else {
				cout << "Error";
			}
		}
		if (point == false) {
			scale = 0;
		}
	}
	else {
		sign = 1;
		bool point = false;
		for (int i = 0; i < numstr.size(); i++) {
			if (numstr[i] == '.' && point == false) {
				point = true;
				scale = numstr.size() - i -1;
				continue;
			}
			else if (isdigit(numstr[i])) {
				int tmp = numstr[i] - '0';
				num.push_back(tmp);
			}
			else {
				cout << "Error";
			}
		}
	}
}

FixedPointNum::FixedPointNum(FixedPointNum & f)
{
	this->num = f.num;
	this->scale = f.scale;
	this->sign = f.sign;
}

FixedPointNum::FixedPointNum(vector<int> n, int s, int ss)
{
	num = n;
	scale = s;
	sign = ss;
}

bool abs_greater_than(FixedPointNum f1, FixedPointNum f2) {
	int scale1 = f1.get_scale();
	int scale2 = f2.get_scale();
	int maxscale = (scale1 > scale2) ? scale1 : scale2;
	f1.change_scale(maxscale);
	f2.change_scale(maxscale);

	vector<int> n1 = f1.get_num();
	vector<int> n2 = f2.get_num();
	int i1 = 0;
	while (n1[i1] == 0 && i1< n1.size()) {
		i1++;
		if (i1 == n1.size()) break;
	}
	int i2 = 0;
	while (n2[i2] == 0 && i2 < n2.size()) {
		i2++;
		if (i2 == n2.size()) break;
	}
	if (n1.size() - i1 > n2.size() - i2) {
		return true;
	}
	else if(n1.size() - i1 < n2.size() - i2){
		return false;
	}
	else {
		while (i1 < n1.size() && i2 < n2.size()) {
			if (n1[i1] > n2[i2]) {
				return true;
			}
			i1++;
			i2++;
		}
		return false;
	}
}

FixedPointNum & FixedPointNum::operator=(FixedPointNum & f)
{
	this->num = f.get_num();
	this->scale = f.get_scale();
	this->sign = f.get_sign();
	return *this;
	// TODO: insert return statement here
}

FixedPointNum & FixedPointNum::operator+=(FixedPointNum & f)
{
	if (this->sign == f.get_sign()) {
		int scale1 = this->get_scale();
		int scale2 = f.get_scale();
		int maxscale = (scale1 > scale2) ? scale1 : scale2;

		this->change_scale(maxscale);
		f.change_scale(maxscale);

		int i1 = this->num.size() - 1;
		int i2 = f.num.size() - 1;
		int carry = 0;

		vector<int> tmpint;
		while (i1 >= 0 && i2 >= 0) {
			int tmp = this->num[i1] + f.num[i2] + carry;
			carry = tmp / 10;
			tmp = tmp % 10;
			tmpint.push_back(tmp);
			i1--;
			i2--;
		}

		while (i1 >= 0) {
			int tmp = this->num[i1] + carry;
			carry = tmp / 10;
			tmp = tmp % 10;
			tmpint.push_back(tmp);
			i1--;
		}
		while (i2 >= 0) {
			int tmp = f.num[i2] + carry;
			carry = tmp / 10;
			tmp = tmp % 10;
			tmpint.push_back(tmp);
			i2--;
		}
		if (carry == 1) {
			tmpint.push_back(1);
		}
		reverse(tmpint.begin(), tmpint.end());
		this->num = tmpint;
		this->scale = maxscale;
		return *this;
	}
	else {
		FixedPointNum fl = (abs_greater_than(*this, f)) ? (*this) : f;
		FixedPointNum fs = (abs_greater_than(*this, f)) ? f : (*this);
		int scale1 = fl.get_scale();
		int scale2 = fs.get_scale();
		int maxscale = (scale1 > scale2) ? scale1 : scale2;
		fl.change_scale(maxscale);
		fs.change_scale(maxscale);
		int totalsign = fl.get_sign();

		vector<int> n1 = fl.get_num();
		vector<int> n2 = fs.get_num();
		int i1 = n1.size() - 1;
		int i2 = n2.size() - 1;
		int carry = 0; // either 0 or -1

		vector<int> tmpint;
		while (i1 >= 0 && i2 >= 0) {
			if (n1[i1]+carry < n2[i2]) {
				int tmp = n1[i1] + carry + 10 - n2[i2];
				carry = -1;
				tmpint.push_back(tmp);
			}
			else {
				int tmp = n1[i1] + carry - n2[i2];
				carry = 0;
				tmpint.push_back(tmp);
			}
			i1--;
			i2--;
		}

		while (i1 >= 0) {
			int tmp = n1[i1] + carry;
			if (tmp >= 0) {
				tmpint.push_back(tmp);
				carry = 0;
			}
			else {
				tmpint.push_back(tmp + 10);
				carry = -1;
			}
			i1--;
		}
		while (i2 >= 0) {
			int tmp = n2[i2] + carry;
			if (tmp >= 0) {
				tmpint.push_back(tmp);
				carry = 0;
			}
			else {
				tmpint.push_back(tmp + 10);
				carry = -1;
			}
			i2--;
		}
		reverse(tmpint.begin(), tmpint.end());
		this->num = tmpint;
		this->scale = maxscale;
		this->sign = totalsign;
		return *this;
	}
// TODO: insert return statement here
}

FixedPointNum FixedPointNum::operator+(FixedPointNum & f)
{
	FixedPointNum result = *this;
	result += f;
	return result;
}

FixedPointNum & FixedPointNum::operator-=(FixedPointNum & f)
{
	FixedPointNum mf(f.get_num(), f.get_scale(), (-1)*f.get_sign());
	*this += mf;
	return *this;
	// TODO: insert return statement here
}

FixedPointNum FixedPointNum::operator-(FixedPointNum & f)
{
	FixedPointNum result = *this;
	result -= f;
	return result;
}

FixedPointNum & FixedPointNum::operator*=(FixedPointNum & f)
{
	FixedPointNum left = *this;
	FixedPointNum sum("0.0");
	int changesign = this->get_sign()*f.get_sign();
	vector<int> n1 = f.get_num();
	for (int i = 0; i < n1.size(); i++) {
		int curr_num = n1[i];
		for (int j = 0; j < curr_num; j++) {
			sum += left;
		}
		if (i != n1.size() - 1)
		{
			sum.num.push_back(0); // similar to *10
		}
	}
	this->num = sum.get_num();
	this->scale = sum.get_scale()+f.get_scale();
	this->sign = changesign;
	return *this;
}

FixedPointNum FixedPointNum::operator*(FixedPointNum & f)
{
	FixedPointNum result = *this;
	result *= f;
	return result;
}

FixedPointNum & FixedPointNum::operator/=(FixedPointNum & f)
{
	FixedPointNum div = *this;
	vector<int> n1 = f.get_num();
	int result_sign = this->sign*f.get_sign();
	int precision = FixedPointNum::get_precision();
	
	int digit_count = 0;
	vector<int> tmpint;
	
	FixedPointNum multiplier("1.0");
	while (abs_greater_than(div,f)) { // end with f>=div
		f.times10();
		multiplier.times10();
	}
	if (!abs_greater_than(div, f) && !abs_greater_than(f, div)) { // div == f
		*this = multiplier;
		this->sign = result_sign;
		return *this;
	}
	// end with f>div

	while (abs_greater_than(f, div)) { // end with div>f but div<10*f
		div.times10();
		multiplier.divide10();
		digit_count++;
		if (digit_count == precision) break;
	}

	// either digit count ==  precision or div>=f

	int result_scale = -1;
	while (digit_count<precision) {
		int curr_num = 0;
		while (abs_greater_than(div, f) ||(!abs_greater_than(div, f) && !abs_greater_than(f,div))) { // div>f , end with div<=f
			div -= f;
			curr_num++;
		}
		tmpint.push_back(curr_num);
		result_scale++;
		div.times10();
		digit_count++;
	}

	if (result_scale == -1) {
		*this = FixedPointNum("0.0");
		return *this;
	}
	else {
		this->num = tmpint;
		this->scale = result_scale;
		this->sign = result_sign;
		*this *= multiplier;
		return *this;
	}
	// TODO: insert return statement here
}

FixedPointNum FixedPointNum::operator/(FixedPointNum & f)
{
	FixedPointNum result = *this;
	result /= f;
	return result;
}

void FixedPointNum::times10()
{
	num.push_back(0);
}

void FixedPointNum::divide10()
{
	vector<int> n = num;
	reverse(n.begin(), n.end());
	n.push_back(0);
	scale++;
	reverse(n.begin(), n.end());
	num = n;
}

void FixedPointNum::print()
{
	if (sign == -1) cout << "-";
	bool zero = true;
	for (int i = 0; i < num.size(); ++i) {
		if (num[i] == 0 && zero == true && (i + scale) <= num.size()-2) continue;
		if (num[i] != 0) zero = false;
		if ((i + scale) == num.size())
			cout << '.';
		cout << num[i];
	}
}

int FixedPointNum::get_scale()
{
	return scale;
}

void FixedPointNum::change_scale(int s)
{
	for (int i = scale; i < s; ++i) {
		num.push_back(0);
	}
	scale = s;
}

int FixedPointNum::get_sign()
{
	return sign;
}

vector<int> FixedPointNum::get_num()
{
	return num;
}


int main()
{

	// 1A
	float f;
	int N;
	string str;
	cout << "1A\nPlease enter a floating point number: ";
	cin >> f;
	cout << "Please enter the number of times you want it to be added:";
	cin >> N;

	float sum = 0;

	for (int i = 0; i<N; i++) {
		sum = sum + f;
	}
	//	 incorrect result for f = 0.0001, N = 10000

	cout << "Result is: " << sum << endl << endl;

	cout << "1B \nUsing my fixed point Code:\n";
	N = 10000;
	int fsum = 0;
	int ff = 1;
	int multiplier = 10000;

	for (int i = 0; i<N; i++) {
		fsum = fsum + ff;
	}
	int intpart = fsum / multiplier;
	int decpart = fsum % multiplier;

	cout << "Result is: " << intpart << ".";
	decpart = decpart * 10;
	while (multiplier != 1) {
		cout << decpart / multiplier;
		decpart = decpart % multiplier;
		multiplier /= 10;
	}
	cout << endl << endl;


	cout << "1C \nSame Problem Using My Fixed Point Number:" << endl;
	FixedPointNum newsum("0.0");
	FixedPointNum b("0.0001");
	
	cout << endl;
	for (int i = 0; i < N; i++) {
		newsum += b;
	}
	cout << "Result is: ";
	newsum.print();
	cout << endl << endl;

	system("pause");

	// test codes

	//FixedPointNum f1("-0.000002");
	//FixedPointNum f2("0.000001");
	//cout << abs_greater_than(f1, f2)<<endl;
	//f1.print();
	//f1 += f2;
	//cout << endl;
	//f1.print();
	//cout << endl;
	//FixedPointNum f3("-1.2");
	//FixedPointNum f5("0.0");
	//f5 += f3;
	//f5.print();
	//f3.print();
	//FixedPointNum f4("1.2");
	//f4.print();
	//f3 *= f3;
	//f3.print();
	//cout << endl;
	//f3 /= f2;
	//f3.print();
	//cout << endl;
	//FixedPointNum f8("1.0");
	//FixedPointNum f9("3.0");
	//f8/= f9;
	//f8.print();
	//
	//// next steps:
	//// operator + overloading
	//// operator >> and << overloading
	//// -,*,/ overloading
	//// error handling
}

