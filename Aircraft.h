#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "utils.h"
#include "SIMPLEX.h"
#include <functional>

class Aircraft {

	// physical parameters
	const double mass = 11; // kg
	const double g = 9.81; // m/s^2
	const double PI = 3.14;  /* pi */
	const double S = 0.55; // wing area, m ^ 2
	const double c_bar = 0.18994; // mean aero chord, m
	const double b = 2.8956; // wing span, m

	// interial values // precalculated 
	const double Ixx = 0.8244;  // kg/m^3
	const double Iyy = 1.135;
	const double Izz = 1.759;
	const double Ixz = 0.1204;
	const double gamma = Ixx * Izz - Ixz * Ixz;
	const double gamma11 = Ixz * (Ixx - Iyy + Izz) / gamma;
	const double gamma12 = (Izz * (Izz - Iyy) + Ixz * Ixz) / gamma;
	const double gamma13 = Izz / gamma;
	const double gamma14 = Ixz / gamma;
	const double gamma21 = (Izz - Ixx) / Izz;
	const double gamma22 = Ixz / Izz;
	const double gamma31 = ((Ixx - Izz) + Ixz * Ixz) / gamma;
	const double gamma32 = Ixx / gamma;

	// {u,v,w,p,q,r,e0,e1,e2,e3,phi,theta,psi,Xe,Ye,Ze}
	std::vector<double> states{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	// {V, alpha, beta}
	std::vector<double> velocities{ 0.0, 0.0, 0.0 };

	// {phi, theta, psi}
	std::vector<double> euler_angles{ 0.0, 0.0, 0.0 };

	// {throttle, elevator, aileron, rudder}
	std::vector<double> controls{ 0.0, 0.0, 0.0, 0.0};
	const double max_elv = 45.0 * PI / 180.0;
	const double max_ail = 30.0 * PI / 180.0;
	const double max_rud = 45.0 * PI / 180.0;

	// Forces // {Fx,Fy,Fz}
	std::vector<double> Forces{ 0.0, 0.0, 0.0 };
	std::vector<double> Forces_aero{ 0.0, 0.0, 0.0 };
	std::vector<double> Forces_prop{ 0.0, 0.0, 0.0 };
	std::vector<double> Forces_grav{ 0.0, 0.0, 0.0 };

	// Moments // {Mx,My,Mz}
	std::vector<double> Moments{ 0.0, 0.0, 0.0 };
	std::vector<double> Moments_aero{ 0.0, 0.0, 0.0 };
	std::vector<double> Moments_prop{ 0.0, 0.0, 0.0 };
	std::vector<double> Moments_grav{ 0.0, 0.0, 0.0 };

	// aerodynamic coeffs.
	double oswald_e = 0.9;
	double alpha_0 = 0.47;
	double M = 50.0; // transition rate// Fixed
	double CL_0 = 0.23;
	double CL_alpha = 5.61;
	double CL_q = 7.95;
	double CL_delta_e = 0.13;
	double CD_0 = 0.043;
	double CD_alpha = 0.03;
	double CD_p = 0.04;
	double CD_q = 0.0;
	double CD_delta_e = 0.0135;
	double CY_0 = 0.0;
	double CY_beta = -0.98;
	double CY_p = 0.0;
	double CY_r = 0.0;
	double CY_delta_a = 0.075;
	double CY_delta_r = 0.19;
	double Cm_0 = 0.0135;
	double Cm_alpha = -2.74;
	double Cm_q = -38.21;
	double Cm_delta_e = -0.99;
	double Cl_0 = 0.0;
	double Cl_beta = -0.13;
	double Cl_p = -0.51;
	double Cl_r = 0.25;
	double Cl_delta_a = 0.17;
	double Cl_delta_r = 0.0024;
	double Cn_0 = 0.0;
	double Cn_beta = 0.073;
	double Cn_p = 0.069;
	double Cn_r = -0.095;
	double Cn_delta_a = -0.011;
	double Cn_delta_r = -0.069;

	// atmospheric paramters
	const double rho0 = 1.2682; // density sea level // kg/m^3
	double rho = 1.2682; // Kg / m3

	// proulsion parameters
	const double D_prop = 20 * (0.0254); // 0.508; // m
	const double Kv = 145.0;
	const double Kq = (1.0 / Kv) * 60.0 / (2.0 * PI); // N*m/A
	const double R_motor = 0.042;  // ohms
	const double i0 = 1.5; // Amp
	const double ncells = 12.0;
	const double Volt_max = 3.7 * ncells;  // max volt
	const double C_Q2 = -0.01664;
	const double C_Q1 = 0.004970;
	const double C_Q0 = 0.005230;
	const double C_T2 = -0.1079;
	const double C_T1 = -0.06044;
	const double C_T0 = 0.09357;

	// writing file param
	bool initFlag = 1;
	double current_time = 0.0;

	// trim
	double trim_velocity = 0.0;
	double trim_altitude = 0.0;
	double trim_gamma = 0.0;
	double trim_radius = 0.0;

public:
	Aircraft() {};
	Aircraft(const std::vector<double>& intial_states);
	~Aircraft() = default;

private:
	void atmosphere_properties();
	double trim_cost_func(std::vector<double> Z);

public:
	void update(const double& time_step);
	void set_control_inputs(const double& throttle, const double& elevator, const double& aileron, const double& rudder);
	void write_to_file(std::ofstream& output_stream, const char* filename = "simulation_result.csv");
	double get_altitude() const;
	std::vector<double> get_states() const;
	std::vector<double> get_controls() const;
	std::vector<double> trim(const double& velocity, const double& flight_path_angle, const double& altitude, const double& orbit_radius);
	std::vector<double> calculate_derivatives(const std::vector<double>& input_states, const std::vector<double>& control_inputs);
	std::vector<double> quat2euler(const double& e0, const double& e1, const double& e2, const double& e3);

private:
	void altitude_check();
	void quat2euler();
	void norm_quaternion();
	void euler2quat(double yaw, double pitch, double roll);
};

