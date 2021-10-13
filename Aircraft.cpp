#include "Aircraft.h"

Aircraft::Aircraft(const std::vector<double>& intial_states) {
	for (int i = 0; i < (int)intial_states.size(); i++) {
		states[i] = intial_states[i];
	}
	
	// check if quaternions are properly initialised
	if (states[6] == 0.0 && abs(sqrt(states[7] * states[7] + states[8] * states[8] + states[9] * states[9]) - 1) > 0.01) {
		states[6] = 1.0;
		states[7] = 0.0;
		states[8] = 0.0;
		states[9] = 0.0;
	}

	quat2euler(); // init euler angles

	velocities[0] = states[0]; // init V
}

void Aircraft::atmosphere_properties() {
	double alt = -states[12] / 3.281;
	if (alt <= 36089.0) {
		rho = rho0 * pow((1 - 0.0000068756 * alt), 5.2559);
	}
	else if (alt <= 65617.0) {
		rho = rho0 * 0.22336 * exp((36089.0 - alt) / 20806);
	}
}

std::vector<double> Aircraft::calculate_derivatives(const std::vector<double>& input_states, const std::vector<double>& control_inputs) {

	double u = input_states[0], v = input_states[1], w = input_states[2];
	double p = input_states[3], q = input_states[4], r = input_states[5];
	double e0 = input_states[6], e1 = input_states[7], e2 = input_states[8], e3 = input_states[9];

	double d_thrtl = control_inputs[0];
	double d_elv = control_inputs[1];
	double d_ail = control_inputs[2];
	double d_rud = control_inputs[3];

	double Va = sqrt(u * u + v * v + w * w);
	double alpha = atan2(w, u);
	double beta = asin(v / Va);

	// calculate aerodynamics forces and moments
	double s_alpha = sin(alpha), c_alpha = cos(alpha);

	// non - dimensional
	double p_nd = 0.5 * p * b / Va;
	double q_nd = 0.5 * q * c_bar / Va;
	double r_nd = 0.5 * r * b / Va;

	double sigmoid_alpha = (1 + exp(-M * (alpha - alpha_0)) + exp(M * (alpha + alpha_0))) / ((1 + exp(-M * (alpha - alpha_0))) * (1 + exp(M * (alpha + alpha_0))));
	double CL_alpha_rec = (1 - sigmoid_alpha) * (CL_0 + CL_alpha * alpha) + sigmoid_alpha * (2 * sign(alpha) * pow(s_alpha, 2) * c_alpha);
	double CD_alpha = CD_p + pow((CL_0 + CL_alpha_rec * alpha), 2) / (PI * oswald_e * b * b / S);
	double CL = CL_alpha_rec + CL_q * q_nd + CL_delta_e * d_elv;
	double CD = CD_alpha + CD_q * q_nd + CD_delta_e * d_elv;

	// aerodynamic forces
	double CX = -c_alpha * CD + s_alpha * CL;
	double CZ = -s_alpha * CD - c_alpha * CL;
	double CY = CY_0 + CY_beta * beta + CY_p * p_nd + CY_r * r_nd + CY_delta_a * d_ail + CY_delta_r * d_rud;

	// aerodynamic moments
	double Cl = Cl_0 + beta * Cl_beta + p_nd * Cl_p + r_nd * Cl_r + d_ail * Cl_delta_a + d_rud * Cl_delta_r;
	double Cm = Cm_0 + alpha * Cm_alpha + q_nd * Cm_q + d_elv * Cm_delta_e;
	double Cn = Cn_0 + beta * Cn_beta + p_nd * Cn_p + r_nd * Cn_r + d_ail * Cn_delta_a + d_rud * Cn_delta_r;

	// dynamic pressure
	double q_bar = 0.5 * rho * Va * Va;

	Forces_aero[0] = q_bar * S * CX;
	Forces_aero[1] = q_bar * S * CY;
	Forces_aero[2] = q_bar * S * CZ;
	Moments_aero[0] = q_bar * S * b * Cl;
	Moments_aero[1] = q_bar * S * c_bar * Cm;
	Moments_aero[2] = q_bar * S * b * Cn;

	// gravity forces
	Forces_grav[0] = mass * g * 2 * (e1 * e3 - e2 * e0);
	Forces_grav[1] = mass * g * 2 * (e2 * e3 + e1 * e0);
	Forces_grav[2] = mass * g * (e3 * e3 + e0 * e0 - e1 * e1 - e2 * e2);

	// propulsion forces
	double a = C_Q0 * (rho * pow(D_prop, 5)) / pow((2 * PI), 2);
	double b = (rho * pow(D_prop, 4) * C_Q1 * Va) / (2 * PI) + (Kq * Kq) / R_motor;
	double c = rho * pow(D_prop, 3) * C_Q2 * (Va * Va) - (Kq * Volt_max * d_thrtl) / R_motor + Kq * i0;

	double omega_op = (-b + sqrt((b * b) - 4 * a * c)) / (2 * a);
	double J_op = (2 * PI * Va) / (omega_op * D_prop);
	double CT = C_T2 * (J_op * J_op) + C_T1 * J_op + C_T0;
	double CQ = C_Q2 * (J_op * J_op) + C_Q1 * J_op + C_Q0;

	Forces_prop[0] = CT * (rho * (omega_op * omega_op) * pow(D_prop, 4)) / pow((2 * PI), 2);
	Moments_prop[0] = rho * pow((omega_op / (2 * PI)), 2) * pow(D_prop, 5) * CQ;

	// Forces and moments 
	Forces[0] = Forces_aero[0] + Forces_grav[0] + Forces_prop[0];
	Forces[1] = Forces_aero[1] + Forces_grav[1];
	Forces[2] = Forces_aero[2] + Forces_grav[2];

	Moments[0] = Moments_aero[0] - Moments_prop[0];
	Moments[1] = Moments_aero[1];
	Moments[2] = Moments_aero[2];

	// calculate derivatives
	std::vector<double> d_states{input_states};
	d_states[0] = r * v - q * w + Forces[0] / mass;
	d_states[1] = p * w - r * u + Forces[1] / mass;
	d_states[2] = q * u - p * v + Forces[2] / mass;
	d_states[3] = gamma11 * p * q - gamma12 * q * r + gamma13 * Moments[0] + gamma14 * Moments[2];
	d_states[4] = gamma21 * p * r - gamma22 * (p * p - r * r) + Moments[1] / Iyy;
	d_states[5] = gamma31 * p * q - gamma11 * q * r + gamma14 * Moments[0] + gamma32 * Moments[2];
	d_states[6] = 0.5 * (-p * e1 - q * e2 - r * e3);
	d_states[7] = 0.5 * (p * e0 + r * e2 - q * e3);
	d_states[8] = 0.5 * (q * e0 - r * e1 + p * e3);
	d_states[9] = 0.5 * (r * e0 + q * e1 - p * e2);
	d_states[10] = u * (e1 * e1 + e0 * e0 - e2 * e2 - e3 * e3) + 2 * v * (e1 * e2 - e3 * e0) + 2 * w * (e1 * e3 + e2 * e0);
	d_states[11] = 2 * u * (e1 * e2 + e3 * e0) + v * (e2 * e2 + e0 * e0 - e1 * e1 - e3 * e3) + 2 * w * (e2 * e3 - e1 * e0);
	d_states[12] = -(-2 * u * (e1 * e3 - e2 * e0) - 2 * v * (e2 * e3 + e1 * e0) - w * (e3 * e3 + e0 * e0 - e1 * e1 - e2 * e2));

	return d_states;
}

void Aircraft::update(const double & time_step) {

	// RK4
	std::vector<double> k1, k2, k3, k4, n_states;
	k1 = calculate_derivatives(states, controls);
	k2 = calculate_derivatives(states + time_step * 0.5 * k1, controls);
	k3 = calculate_derivatives(states + time_step * 0.5 * k2, controls);
	k4 = calculate_derivatives(states + time_step * k3, controls);
	n_states = (k1 + 2 * k2 + 2 * k3 + k4) * (time_step / 6.0);

	states[0] += n_states[0];
	states[1] += n_states[1];
	states[2] += n_states[2];
	states[3] += n_states[3];
	states[4] += n_states[4];
	states[5] += n_states[5];
	states[6] += n_states[6];
	states[7] += n_states[7];
	states[8] += n_states[8];
	states[9] += n_states[9];
	norm_quaternion();
	quat2euler();
	states[10] += n_states[10];
	states[11] += n_states[11];
	states[12] += n_states[12];

	current_time += time_step;
	// update velocities
	velocities[0] = sqrt(states[0] * states[0] + states[1] * states[1] + states[2] * states[2]);// m / sec
	states[0] ? (velocities[1] = atan2(states[2], states[0])) : (velocities[1] = sign(states[0]) * PI / 2.0);// rad
	velocities[0] ? (velocities[2] = asin(states[1] / velocities[0])) : (velocities[2] = sign(velocities[0]) * PI / 2.0);// rad
}

void Aircraft::set_control_inputs(const double& throttle, const double& elevator, const double& aileron, const double& rudder) {
	// inputs in rad // not in degrees
	if (throttle < 0.0) controls[0] = 0.0;
	else if (throttle > 1.0) controls[0] = 1.0;
	else controls[0] = throttle;

	if (abs(elevator) <= 1e-3) controls[1] = 0.0;
	else if (elevator < -max_elv) controls[1] = -max_elv;
	else if (elevator > max_elv) controls[1] = max_elv;
	else controls[1] = elevator;// *PI / 180.0;

	if (abs(aileron) <= 1e-3) controls[2] = 0.0;
	else if (aileron < -max_ail) controls[2] = -max_ail;
	else if (aileron > max_ail) controls[2] = max_ail;
	else controls[2] = aileron;// *PI / 180.0;

	if (abs(rudder) <= 1e-3) controls[3] = 0.0;
	else if (rudder < -max_rud) controls[3] = -max_rud;
	else if (rudder > max_rud) controls[3] = max_rud;
	else controls[3] = rudder;// *PI / 180.0;
}

void Aircraft::altitude_check() { // Terminate the simulation
	if ((-states[12]) < 0.0 || (-states[12]) > 20000.0) {
		std::cout << "Altitude is below 0 Km or above 20 km" << std::endl;
		std::exit(EXIT_SUCCESS);
	}
};

void Aircraft::quat2euler() {
	const double e0 = states[6], e1 = states[7], e2 = states[8], e3 = states[9];
	euler_angles[0] = atan2(2 * (e0 * e1 + e2 * e3), e0 * e0 + e3 * e3 - e1 * e1 - e2 * e2);
	euler_angles[1] = asin(2 * (e0 * e2 - e1 * e3));
	euler_angles[2] = atan2(2 * (e0 * e3 + e1 * e2), e0 * e0 + e1 * e1 - e2 * e2 - e3 * e3);
}

std::vector<double> Aircraft::quat2euler(const double& e0, const double& e1, const double& e2, const double& e3) {
	std::vector<double> euler_angles{ 0.0, 0.0, 0.0 };
	euler_angles[0] = atan2(2 * (e0 * e1 + e2 * e3), e0 * e0 + e3 * e3 - e1 * e1 - e2 * e2);
	euler_angles[1] = asin(2 * (e0 * e2 - e1 * e3));
	euler_angles[2] = atan2(2 * (e0 * e3 + e1 * e2), e0 * e0 + e1 * e1 - e2 * e2 - e3 * e3);
	return euler_angles;
}

void Aircraft::norm_quaternion() {
	double normE = sqrt(states[6] * states[6] + states[7] * states[7] + states[8] * states[8] + states[9] * states[9]);
	states[6] /= normE;
	states[7] /= normE;
	states[8] /= normE;
	states[9] /= normE;
}

void Aircraft::write_to_file(std::ofstream& output_stream, const char* filename) {

	if (initFlag) {
		output_stream.open(filename, std::ios_base::out | std::ios_base::trunc);
		output_stream << "t" 
		<< "," << "u" << "," << "v" << "," << "w" << ","
		<< "p" << "," << "q" << "," << "r" << ","
		<< "e0" << "," << "e1" << "," << "e2" << "," << "e3" << ","
		<< "Xe" << "," << "Ye" << "," << "h" << ","
		<< "throttle" << "," << "elevator" << ","
		<< "aileron" << "," << "rudder" << ","
		<< "Fx" << "," << "Fy" << "," << "Fz" << ","
		<< "Mx" << "," << "My" << "," << "Mz" << std::endl;
		initFlag = 0;
	}

	output_stream << current_time << "," << states[0] << "," << states[1] << "," << states[2] << ","
		<< states[3] << "," << states[4] << "," << states[5] << ","
		<< states[6] << "," << states[7] << "," << states[8] << "," << states[9] << ","
		<< states[10] << "," << states[11] << "," << (-states[12]) << ","
		<< controls[0] << "," << controls[1] << ","
		<< controls[2] << "," << controls[3] << ","
		<< Forces[0] << "," << Forces[1] << "," << Forces[2] << ","
		<< Moments[0] << "," << Moments[1] << "," << Moments[2] << std::endl;
}

double Aircraft::get_altitude() const {
	return -states[12];
}

std::vector<double> Aircraft::get_states() const {
	return states;
}

std::vector<double> Aircraft::get_controls() const {
	return controls;
}

void Aircraft::euler2quat(double yaw, double pitch, double roll) // yaw (Z), pitch (Y), roll (X)
{
	// Abbreviations for the various angular functions
	double cy = cos(yaw * 0.5);
	double sy = sin(yaw * 0.5);
	double cp = cos(pitch * 0.5);
	double sp = sin(pitch * 0.5);
	double cr = cos(roll * 0.5);
	double sr = sin(roll * 0.5);
	
	states[6] = cr * cp * cy + sr * sp * sy;
	states[7] = sr * cp * cy - cr * sp * sy;
	states[8] = cr * sp * cy + sr * cp * sy;
	states[9] = cr * cp * sy - sr * sp * cy;
}

std::vector<double> Aircraft::trim(const double& velocity, const double& flight_path_angle, const double& altitude, const double& orbit_radius) {

	// initial states
	trim_gamma = flight_path_angle;
	trim_velocity = velocity;
	trim_altitude = altitude;
	trim_radius = orbit_radius;
	euler2quat(0.0, trim_gamma, 0.0); // yaw, pitch, roll

	std::vector<double> Z; // 13 states and 4 controls
	for (int i = 0; i < 13; i++) { Z.push_back(states[i]); }
	for (int i = 0; i < 4; i++) { Z.push_back(controls[i]); }

	std::vector<double> Zstar;
	double temp_va, temp_h;

	int iter = 0;
	const int max_try = 10; // 10;
	while (true) {
		if (iter >= max_try) break;
		if (trim_cost_func(Z) <= 1e-10) break;
		Zstar = SIMPLEX(std::bind(&Aircraft::trim_cost_func, this, std::placeholders::_1), Z, 1e-10, 1e-10, 10000);
		temp_va = sqrt(Zstar[0] * Zstar[0] + Zstar[1] * Zstar[1] + Zstar[2] * Zstar[2]);
		temp_h = -Zstar[12];
		if (abs(temp_va - trim_altitude) <= 1e-1 && abs(temp_h - trim_altitude) <= 1e-1) {
			break;
		}
		iter += 1;
		Z = Zstar;
	}

	for (int i = 0; i < 3; i++) { states[i] = abs(Zstar[i]) >= 1e-2 ? Zstar[i] : 0.0; }
	for (int i = 3; i < 10; i++) { states[i] = abs(Zstar[i]) >= 1e-3 ? Zstar[i] : 0.0; }
	states[10] = 0.0; // Xe
	states[11] = 0.0;  // Ye
	states[12] = -altitude; // trim_altitude; // Zstar[12]; // -he
	norm_quaternion();
	quat2euler();
	set_control_inputs(Zstar[13], Zstar[14], Zstar[15], Zstar[16]);
	return Zstar;
}

double Aircraft::trim_cost_func(std::vector<double> Z) {

	std::vector<double> X{ Z.begin(), Z.begin() + 13 };
	std::vector<double> U{ Z.begin() + 13, Z.end() };

	std::vector<double> Xdot = calculate_derivatives(X, U);

	double Va = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
	std::vector<double> eulers = quat2euler(X[6], X[7], X[8], X[9]);

	std::vector<double> constraints;
	constraints.push_back(Xdot[0]); // udot
	constraints.push_back(Xdot[1]); // vdot
	constraints.push_back(Xdot[2]); // wdot
	constraints.push_back(Xdot[3]); // pdot
	constraints.push_back(Xdot[4]); // qdot
	constraints.push_back(Xdot[5]); // rdot
	constraints.push_back(Xdot[7]); // e1dot
	constraints.push_back(Xdot[8]); // e2dot
	constraints.push_back(Va - trim_velocity);
	constraints.push_back(X[1]); // v
	constraints.push_back(X[6] * X[6] + X[7] * X[7] + X[8] * X[8] + X[9] * X[9] - 1.0); // quat = 1;
	constraints.push_back(-Va * sin(trim_gamma) - Xdot[12]); // gamma effect on altitude

	// only add if trim is for steady level or climbing flight
	if (trim_radius == 0.0) {
		constraints.push_back(Xdot[6]); // e0dot
		constraints.push_back(Xdot[9]); // e3dot // for psi-dot
		constraints.push_back(X[3]); // p
		constraints.push_back(X[4]); // q
		constraints.push_back(X[5]); // r
		constraints.push_back(X[7]); // e1 => phi
		constraints.push_back(X[9]); // e3 => psi
		constraints.push_back(X[12] + trim_altitude);
	}
	else {
		double psi_dot = X[4] * sin(eulers[0]) / cos(eulers[1]) + X[5] * cos(eulers[0])  / cos(eulers[1]);
		constraints.push_back(psi_dot - cos(trim_gamma) * trim_velocity / trim_radius); // psi-dot = Va/R
	}

	double cost = 0.0;
	for (double& i : constraints) cost += (i * i);

	return cost;
}
