#pragma once

/**
* Use:
*  PID myPID;  // init pid
*  myPID.set_gains(kp,ki,kd);
*  myPID.set_input_range(up_limit, low_limit); // saturation limit of the input
*  input = myPID.caclulate(y_desired, y_actual);   // calculate inputs required
*          "or"
*  input = myPID.calculate_digital(y_desired, y_actual); // digital implementation of PID
*
*/
class PID {
public:
	PID();

	virtual ~PID();

	/**
	* set kp, ki and kd
	* @param kp Proportional Gain
	* @param ki Integral Gain
	* @param kd Derivative Gain
	*/
	void set_gains(double kp, double ki, double kd);

	/**
	* set input limits (saturation)
	* @param up_limit Upper Saturation limit
	* @param low_limit lower Saturation limit
	*/
	void set_input_range(double up_limit, double low_limit);

	/**
	* simple pid contoller
	* @param y_d Desired output
	* @param y Current Output
	* @return Output Command
	*/
	double calculate(double y_d, double y);

	/**
	* Digital implementation of PID
	* @param y_d Desired output
	* @param y Current Output
	* @return Output Command
	* @note Equations for digital implementation of PID can be found in below reference.
	* @cite {"Digital Implementation of PID Controller for Temperature Control", PrachiRusia,
	*      International Journal of Scientific & Engineering Research Volume 8, Issue 5, May-2017,
	*      ISSN 2229-5518}
	*/
	double calculate_digital(double y_d, double y);

private:
	// for both implementation
	double output{ 0 };
	double kp{ 1 };
	double ki{ 0 };
	double kd{ 0 };
	double up_limit{ 1 };
	double low_limit{ 0 };

	// for simple pid
	double I{ 0 }; // integrator
	double I_max{ 0 };    // for anti-windup
	double D; // differentiator
	double prev_error;
	double dt = 0.01;   // time step    // hard coded for faster code

	// for digital pid
	double error_d1; // for digital implementation
	double error_d2;
	double y_d1{ 0 };
	double K1{ 0 };
	double K2{ 0 };

	/**
	* check for saturation of the output
	* @param output unsaturated output
	* @return output saturated output
	*/
	double saturation_check(double output);

	/**
	* check anti-windup for integrator
	* @param error error between desired and actual value
	* @return error_ appropriate error for integral only
	*/
	double integrator_anti_windup(double error);
};
