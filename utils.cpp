#include "utils.h"

std::vector<double> operator+(const std::vector<double> &lhs, const std::vector<double> &rhs)
{
	if (lhs.size() != rhs.size())
	{ // size check
		throw std::runtime_error("Can't add two vectors of different sizes");
	}

	std::vector<double> result;
	for (int i = 0; i < (int)lhs.size(); i++)
	{
		result.push_back(lhs.at(i) + rhs.at(i));
	}
	return result;
}

std::vector<double> operator-(const std::vector<double> &lhs, const std::vector<double> &rhs)
{
	if (lhs.size() != rhs.size())
	{ // size check
		throw std::runtime_error("Can't add two vectors of different sizes");
	}

	std::vector<double> result;
	for (int i = 0; i < (int)lhs.size(); i++)
	{
		result.push_back(lhs.at(i) - rhs.at(i));
	}
	return result;
}

std::vector<double> operator*(const double &n, const std::vector<double> &rhs)
{
	std::vector<double> result;
	for (int i = 0; i < (int)rhs.size(); i++)
	{
		result.push_back(rhs.at(i) * n);
	}
	return result;
}

std::vector<double> operator*(const std::vector<double> &lhs, const double &n)
{
	std::vector<double> result;
	for (int i = 0; i < (int)lhs.size(); i++)
	{
		result.push_back(lhs.at(i) * n);
	}
	return result;
}

std::vector<double> operator/(const std::vector<double> &lhs, const double &n)
{
	std::vector<double> result;
	for (int i = 0; i < (int)lhs.size(); i++)
	{
		result.push_back(lhs.at(i) / n);
	}
	return result;
}

std::vector<double> abs(const std::vector<double> &v)
{
	std::vector<double> other = v;
	for (double &i : other)
	{
		i >= 0 ? other.push_back(i) : other.push_back(-i);
	}
	return other;
}

double max(const std::vector<double> &v)
{
	return *std::max_element(std::begin(v), std::end(v));
}

double sign(double &a)
{
	if (a >= 0.0)
	{
		return 1.0;
	}
	else
	{
		return -1.0;
	}
}
