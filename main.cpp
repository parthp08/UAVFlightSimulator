#include <iostream>
#include <fstream>
#include "Aircraft.h"
#include "PID.h"

constexpr double PI = 3.14;
constexpr double DEGREE2RAD = 0.01745;

void print_vector(std::vector<double>& vec) {
	for (double& i : vec) std::cout << "\t" << i;
	std::cout << std::endl;
}


int main() {

	std::cout << "\t WELCOME TO FLIGHT SIMULATION \t\n";

	// inital states
	// states =>			      u,    v,   w,  p,   q,   r,    e0, e1,  e2,  e3,  Xe,  Ye,   Ze
	std::vector<double> states { 17.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -100.0 };

	// intialise aircraft from intial states
	Aircraft uav(states);
	
	// simulation parameters
	double t_current = 0.0;
	double t_end = 60.0;
	double t_step = 1e-2;

	//// altitude hold controller setup
	//PID altitude_hold_controller;
	//const double Kp = 0.2;
	//const double Ki = 0.05;
	//const double Kd = 0.06;
	//altitude_hold_controller.set_gains(Kp, Ki, Kd);
	//altitude_hold_controller.set_input_range(45 * PI / 180.0, -45 * PI / 180.0);

		// intial controls
	double throttle = 0.5, evelvator = 0.0, aileron = 0.0, rudder = 0.0;
	uav.set_control_inputs(throttle, evelvator, aileron, rudder);

	const double desired_velocity = 17.0;
	const double desired_flight_angle = 0.0 * DEGREE2RAD; // rad
	const double desired_altitude = 100.0;
	const double desired_turning_radius = 100.0; // 0.0 => no turn
	std::vector<double> trim = uav.trim(desired_velocity, desired_flight_angle, desired_altitude, desired_turning_radius);

	// trim states and controls
	std::vector<double> Xstar, Ustar;
	Xstar = uav.get_states();
	Ustar = uav.get_controls();

	std::cout << "\ttrim states : ";
	for (auto i : Xstar) { std::cout << "\t" << i; }
	std::cout << std::endl;

	std::cout << "\ttrim controls : ";
	for (auto i : Ustar) { std::cout << "\t" << i; }
	std::cout << std::endl;

	std::cout << " - Initiating Flight Simulation" << std::endl;

	std::cout << "verification of trim : " << std::endl;
	std::vector<double> XdotStar = uav.calculate_derivatives(Xstar, Ustar);
	double VaStar = sqrt(Xstar[0] * Xstar[0] + Xstar[1] * Xstar[1] + Xstar[2] * Xstar[2]);
	double vStar = Xstar[1];
	std::vector<double> eulerStar = uav.quat2euler(Xstar[6], Xstar[7], Xstar[8], Xstar[9]);
	double phiStar = eulerStar[0];
	double thetaStar = eulerStar[1];
	double psiStar = eulerStar[2];
	double gamStar = thetaStar - (atan2(Xstar[2], Xstar[0]));

	std::cout << "XdotStar = ";
	print_vector(XdotStar);
	std::cout << "VaStar = " << VaStar << std::endl;
	std::cout << "gamStar = " << gamStar << std::endl;
	std::cout << "vStar = " << vStar << std::endl;
	std::cout << "phiStar = " << phiStar << std::endl;
	std::cout << "thetaStar = " << thetaStar << std::endl;
	std::cout << "psiStar = " << psiStar << std::endl;

	//// saved data to save time
	////  --- NEED SET STATE FUNCTION TO APPLY THESE ----------------
	//std::vector<double> Xstar{ 24.969, 0.0, 1.24376, 0.0, 0.0, 0.0, 0.99969, 0.0, 0.0248829, 0.0, 39.2636, -0.332861, -100 };
	//std::vector<double> Ustar{ 0.757636, -0.124113, 0.0, 0.0 };

	 //Write data to file stream
	std::ofstream output_stream;

	while (t_current <= t_end && uav.get_altitude() > 0.0) {
		
		uav.update(t_step);

		//elevator = altitude_hold_controller.calculate(-desired_altitude, -uav.get_altitude());
		//uav.set_control_inputs(throttle, elevator, aileron, rudder);

		//if (t_current >= 10 && t_current <= 11) {
		//	uav.set_control_inputs(Ustar[0], Ustar[1]+0.2, Ustar[2], Ustar[3]);
		//}
		//else if (t_current > 11 && t_current <= 12) {
		//	uav.set_control_inputs(Ustar[0], Ustar[1]-0.2, Ustar[2], Ustar[3]);
		//}
		//else {
		//	uav.set_control_inputs(Ustar[0], Ustar[1], Ustar[2], Ustar[3]);
		//}

		uav.write_to_file(output_stream);

		t_current += t_step;
	};

	output_stream.close();
	std::cout << "Flight simultaion complete, please run 'plots.py' in separate console to see the outputs.\n";

	return 0;
}
