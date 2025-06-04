/****************************************************************************
 *
 *   Copyright (c) 2018 - 2019 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file PositionControl.cpp
 */

#include "PositionControl.hpp"
#include "ControlMath.hpp"
#include <float.h>
#include <mathlib/mathlib.h>
#include <px4_platform_common/defines.h>
#include <geo/geo.h>

using namespace matrix;

const trajectory_setpoint_s PositionControl::empty_trajectory_setpoint = {0, {NAN, NAN, NAN}, {NAN, NAN, NAN}, {NAN, NAN, NAN}, {NAN, NAN, NAN}, NAN, NAN};

void PositionControl::setVelocityGains(const Vector3f &P, const Vector3f &I, const Vector3f &D)
{
	_gain_vel_p = P;
	_gain_vel_i = I;
	_gain_vel_d = D;
}

void PositionControl::setVelocityLimits(const float vel_horizontal, const float vel_up, const float vel_down)
{
	_lim_vel_horizontal = vel_horizontal;
	_lim_vel_up = vel_up;
	_lim_vel_down = vel_down;
}

void PositionControl::setThrustLimits(const float min, const float max)
{
	// make sure there's always enough thrust vector length to infer the attitude
	_lim_thr_min = math::max(min, 10e-4f);
	_lim_thr_max = max;
}

void PositionControl::setHorizontalThrustMargin(const float margin)
{
	_lim_thr_xy_margin = margin;
}

void PositionControl::updateHoverThrust(const float hover_thrust_new)
{
	// Given that the equation for thrust is T = a_sp * Th / g - Th
	// with a_sp = desired acceleration, Th = hover thrust and g = gravity constant,
	// we want to find the acceleration that needs to be added to the integrator in order obtain
	// the same thrust after replacing the current hover thrust by the new one.
	// T' = T => a_sp' * Th' / g - Th' = a_sp * Th / g - Th
	// so a_sp' = (a_sp - g) * Th / Th' + g
	// we can then add a_sp' - a_sp to the current integrator to absorb the effect of changing Th by Th'
	const float previous_hover_thrust = _hover_thrust;
	setHoverThrust(hover_thrust_new);

	_vel_int(2) += (_acc_sp(2) - CONSTANTS_ONE_G) * previous_hover_thrust / _hover_thrust
		       + CONSTANTS_ONE_G - _acc_sp(2);
}

void PositionControl::setState(const PositionControlStates &states)
{
	_pos = states.position;
	_vel = states.velocity;
	_yaw = states.yaw;
	_vel_dot = states.acceleration;
}

void PositionControl::setInputSetpoint(const trajectory_setpoint_s &setpoint)
{
	_pos_sp = Vector3f(setpoint.position);
	_vel_sp = Vector3f(setpoint.velocity);
	_acc_sp = Vector3f(setpoint.acceleration);
	_yaw_sp = setpoint.yaw;
	_yawspeed_sp = setpoint.yawspeed;
}

bool PositionControl::update(const float dt)
{
	bool valid = _inputValid();

	if (valid) {
		// _positionControl();
		// _velocityControl(dt);
		SMC_control(dt);

		_yawspeed_sp = PX4_ISFINITE(_yawspeed_sp) ? _yawspeed_sp : 0.f;
		_yaw_sp = PX4_ISFINITE(_yaw_sp) ? _yaw_sp : _yaw; // TODO: better way to disable yaw control
	}

	// There has to be a valid output acceleration and thrust setpoint otherwise something went wrong
	return valid && _acc_sp.isAllFinite() && _thr_sp.isAllFinite();
}

void PositionControl::SMC_control(const float dt)
{
	// Gain of Controller
	float 	k1x = 2,k1y = 2,k1z = 5.55;
	float	k2x = 2,k2y = 2,k2z = 3.2;
	float	c1x = 1.5, c1y=1.5, c1z = 5;

	// Gain of Observer
	float 	lamda1x = 6,lamda2x = 20,lamda3x = 5;
	float	lamda1y = 6,lamda2y = 20,lamda3y = 5;

	// compute error
	// sigma_0 ... estimation error
	ControlMath::setZeroIfNanVector3f(xp_x);
	ControlMath::setZeroIfNanVector3f(xp_y);

	if (!PX4_ISFINITE(Uxx) || !PX4_ISFINITE(Uyy) ){
			Uxx = 0.0f;
			Uyy = 0.0f;
		}

	float U_comp_x = Uxx ;
	float U_comp_y = Uyy ;

	float sigma_0_x = _pos(0) - xp_x(0);
	float sigma_0_y = _pos(1) - xp_y(0);

	float sat_sigmax = ControlMath::sat(sigma_0_x/0.5f);
	float sat_sigmay = ControlMath::sat(sigma_0_y/0.5f);

	xp_x(0) = xp_x(0) + dt*lamda1x *cbrtf(fabsf(sigma_0_x)*fabsf(sigma_0_x))*sign(sigma_0_x) + (dt*xp_x(1) + (dt*dt/2)*xp_x(2)) ;
	xp_x(1) = xp_x(1) + dt*lamda2x*cbrtf(fabsf(sigma_0_x))*sign(sigma_0_x) + (dt*xp_x(2)) +  dt*(U_comp_x);
	xp_x(2) = xp_x(2) + dt*lamda3x*sat_sigmax;

	xp_y(0) = xp_y(0) + dt*lamda1y *cbrtf(fabsf(sigma_0_y)*fabsf(sigma_0_y)*1.0f)*sign(sigma_0_y) + (dt*xp_y(1) + (dt*dt/2)*xp_y(2)) ;
	xp_y(1) = xp_y(1) + dt*lamda2y*cbrtf(fabsf(sigma_0_y))*sign(sigma_0_y) + (dt*xp_y(2)) +  dt*(U_comp_y);
	xp_y(2) = xp_y(2) + dt*lamda3y*sat_sigmay;

	float error_x =_pos(0) -_pos_sp(0);
	float error_y =_pos(1) -_pos_sp(1);
	float error_z =_pos(2) -_pos_sp(2);

	float sx= c1x*error_x +(xp_x(1) -_vel_sp(0)) ;
	float sy= c1y*error_y +(xp_y(1) -_vel_sp(1)) ;
	// z use only Super-Twisting SMC
	float sz= c1z*error_z +(_vel(2) -_vel_sp(2)) ;

	float sat_sx = ControlMath::sat(sx);
	float sat_sy = ControlMath::sat(sy);
	float sat_sz = ControlMath::sat(sz);

	float sat_ex = ControlMath::sat((_pos(0) - xp_x(0)));
	float sat_ey = ControlMath::sat((_pos(1) - xp_y(0)));
	// float sat_ez = ControlMath::sat(_pos(2) -_pos_sp(2));

	sat_sx_int += k2x*sat_sx*dt;
	sat_sy_int += k2y*sat_sy*dt;
	sat_sz_int += k2z*sat_sz*dt;

	// sat_ex_int += sat_ex*dt;
	// sat_ey_int += sat_ey*dt;
	// sat_ez_int += sat_ez*dt;

	float Ux = 0;
	float Uy = 0;
	float Uz = 0;


	// Ux =  c1x*(-xp_x(1) +_vel_sp(0)) - xp_x(2) - lamda2x*cbrt(abs(xp_x(0) -_pos_sp(0)))*sat_ex -k1x*sqrt(abs(sx))*sat_sx - sat_sx_int ;
	// Uy =  c1y*(-xp_y(1) +_vel_sp(1)) - xp_y(2) - lamda2y*cbrt(abs(xp_y(0) -_pos_sp(1)))*sat_ey -k1y*sqrt(abs(sy))*sat_sy - sat_sy_int;
	Ux =  c1x*(-xp_x(1) +_vel_sp(0)) - xp_x(2) - lamda2x*cbrtf(fabsf(xp_x(0) -_pos_sp(0)))*sat_ex -k1x*sqrtf(fabsf(sx))*sat_sx - sat_sx_int ;
	Uy =  c1y*(-xp_y(1) +_vel_sp(1)) - xp_y(2) - lamda2y*cbrtf(fabsf(xp_y(0) -_pos_sp(1)))*sat_ey -k1y*sqrtf(fabsf(sy))*sat_sy - sat_sy_int;
	Uz = c1z*(-_vel(2) +_vel_sp(2)) - k1z*sqrtf(fabsf(sz))*sat_sz - sat_sz_int;

	Uxx = Ux;
	Uyy = Uy;

		// All nan if not takeoff or have a mission
	Vector3f U = Vector3f(Ux, Uy, Uz);
	ControlMath::addIfNotNanVector3f(_acc_sp, U);

	/* Only for Simulation */
	// FILE *fichier = fopen("Obser.txt","a");
	// fprintf(fichier,"  %f\t  %f\t  %f\t   %f\t  %f\t   %f\t    %f\t   %f\t 	%f\t  %f\t   %f\t    %f\t   %f\t   %f\t  %f\t   %f\t    %f\t   %f\t   %f\n",(double)xp_x(1) , (double)_vel(0), (double)xp_y(1),
	// 							(double)_vel(1) , (double)_vel(2) ,(double)xp_x(2),(double)xp_y(2),(double)Ux,
	// 							(double)Uy,(double)Uz,(double)_pos_sp(0), (double)_pos_sp(1), (double)_pos_sp(2),
	// 							(double)_pos(0), (double)_pos(1), (double)_pos(2),
	// 							(double)_vel_sp(0), (double)_vel_sp(1), (double)_vel_sp(2));

	// fclose(fichier);
	/* ------------------- */

	// Call acceleration control to convert to thrust
	_accelerationControl();

	// Saturate maximal vertical thrust
	_thr_sp(2) = math::max(_thr_sp(2), -_lim_thr_max);

	// Get allowed horizontal thrust after prioritizing vertical control
	const float thrust_max_squared = _lim_thr_max * _lim_thr_max;
	const float thrust_z_squared = _thr_sp(2) * _thr_sp(2);
	const float thrust_max_xy_squared = thrust_max_squared - thrust_z_squared;
	float thrust_max_xy = 0;

	if (thrust_max_xy_squared > 0) {
		thrust_max_xy = sqrtf(thrust_max_xy_squared);
	}

	// Saturate thrust in horizontal direction
	const Vector2f thrust_sp_xy(_thr_sp);
	const float thrust_sp_xy_norm = thrust_sp_xy.norm();


	if (thrust_sp_xy_norm > thrust_max_xy) {
		_thr_sp.xy() = thrust_sp_xy / thrust_sp_xy_norm * thrust_max_xy;
	}
}

void PositionControl::_positionControl()
{
	// P-position controller
	Vector3f vel_sp_position = (_pos_sp - _pos).emult(_gain_pos_p);
	// Position and feed-forward velocity setpoints or position states being NAN results in them not having an influence
	ControlMath::addIfNotNanVector3f(_vel_sp, vel_sp_position);
	// make sure there are no NAN elements for further reference while constraining
	ControlMath::setZeroIfNanVector3f(vel_sp_position);

	// Constrain horizontal velocity by prioritizing the velocity component along the
	// the desired position setpoint over the feed-forward term.
	_vel_sp.xy() = ControlMath::constrainXY(vel_sp_position.xy(), (_vel_sp - vel_sp_position).xy(), _lim_vel_horizontal);
	// Constrain velocity in z-direction.
	_vel_sp(2) = math::constrain(_vel_sp(2), -_lim_vel_up, _lim_vel_down);
}

void PositionControl::_velocityControl(const float dt)
{
	// Constrain vertical velocity integral
	_vel_int(2) = math::constrain(_vel_int(2), -CONSTANTS_ONE_G, CONSTANTS_ONE_G);

	// PID velocity control
	Vector3f vel_error = _vel_sp - _vel;
	Vector3f acc_sp_velocity = vel_error.emult(_gain_vel_p) + _vel_int - _vel_dot.emult(_gain_vel_d);

	// No control input from setpoints or corresponding states which are NAN
	ControlMath::addIfNotNanVector3f(_acc_sp, acc_sp_velocity);

	_accelerationControl();

	// Integrator anti-windup in vertical direction
	if ((_thr_sp(2) >= -_lim_thr_min && vel_error(2) >= 0.f) ||
	    (_thr_sp(2) <= -_lim_thr_max && vel_error(2) <= 0.f)) {
		vel_error(2) = 0.f;
	}

	// Prioritize vertical control while keeping a horizontal margin
	const Vector2f thrust_sp_xy(_thr_sp);
	const float thrust_sp_xy_norm = thrust_sp_xy.norm();
	const float thrust_max_squared = math::sq(_lim_thr_max);

	// Determine how much vertical thrust is left keeping horizontal margin
	const float allocated_horizontal_thrust = math::min(thrust_sp_xy_norm, _lim_thr_xy_margin);
	const float thrust_z_max_squared = thrust_max_squared - math::sq(allocated_horizontal_thrust);

	// Saturate maximal vertical thrust
	_thr_sp(2) = math::max(_thr_sp(2), -sqrtf(thrust_z_max_squared));

	// Determine how much horizontal thrust is left after prioritizing vertical control
	const float thrust_max_xy_squared = thrust_max_squared - math::sq(_thr_sp(2));
	float thrust_max_xy = 0.f;

	if (thrust_max_xy_squared > 0.f) {
		thrust_max_xy = sqrtf(thrust_max_xy_squared);
	}

	// Saturate thrust in horizontal direction
	if (thrust_sp_xy_norm > thrust_max_xy) {
		_thr_sp.xy() = thrust_sp_xy / thrust_sp_xy_norm * thrust_max_xy;
	}

	// Use tracking Anti-Windup for horizontal direction: during saturation, the integrator is used to unsaturate the output
	// see Anti-Reset Windup for PID controllers, L.Rundqwist, 1990
	const Vector2f acc_sp_xy_produced = Vector2f(_thr_sp) * (CONSTANTS_ONE_G / _hover_thrust);
	const float arw_gain = 2.f / _gain_vel_p(0);

	// The produced acceleration can be greater or smaller than the desired acceleration due to the saturations and the actual vertical thrust (computed independently).
	// The ARW loop needs to run if the signal is saturated only.
	const Vector2f acc_sp_xy = _acc_sp.xy();
	const Vector2f acc_limited_xy = (acc_sp_xy.norm_squared() > acc_sp_xy_produced.norm_squared())
					? acc_sp_xy_produced
					: acc_sp_xy;
	vel_error.xy() = Vector2f(vel_error) - arw_gain * (acc_sp_xy - acc_limited_xy);

	// Make sure integral doesn't get NAN
	ControlMath::setZeroIfNanVector3f(vel_error);
	// Update integral part of velocity control
	_vel_int += vel_error.emult(_gain_vel_i) * dt;
}

void PositionControl::_accelerationControl()
{
	// Assume standard acceleration due to gravity in vertical direction for attitude generation
	float z_specific_force = -CONSTANTS_ONE_G;

	if (!_decouple_horizontal_and_vertical_acceleration) {
		// Include vertical acceleration setpoint for better horizontal acceleration tracking
		z_specific_force += _acc_sp(2);
	}

	Vector3f body_z = Vector3f(-_acc_sp(0), -_acc_sp(1), -z_specific_force).normalized();
	ControlMath::limitTilt(body_z, Vector3f(0, 0, 1), _lim_tilt);
	// Convert to thrust assuming hover thrust produces standard gravity
	const float thrust_ned_z = _acc_sp(2) * (_hover_thrust / CONSTANTS_ONE_G) - _hover_thrust;
	// Project thrust to planned body attitude
	const float cos_ned_body = (Vector3f(0, 0, 1).dot(body_z));
	const float collective_thrust = math::min(thrust_ned_z / cos_ned_body, -_lim_thr_min);
	_thr_sp = body_z * collective_thrust;
}

bool PositionControl::_inputValid()
{
	bool valid = true;

	// Every axis x, y, z needs to have some setpoint
	for (int i = 0; i <= 2; i++) {
		valid = valid && (PX4_ISFINITE(_pos_sp(i)) || PX4_ISFINITE(_vel_sp(i)) || PX4_ISFINITE(_acc_sp(i)));
	}

	// x and y input setpoints always have to come in pairs
	valid = valid && (PX4_ISFINITE(_pos_sp(0)) == PX4_ISFINITE(_pos_sp(1)));
	valid = valid && (PX4_ISFINITE(_vel_sp(0)) == PX4_ISFINITE(_vel_sp(1)));
	valid = valid && (PX4_ISFINITE(_acc_sp(0)) == PX4_ISFINITE(_acc_sp(1)));

	// For each controlled state the estimate has to be valid
	for (int i = 0; i <= 2; i++) {
		if (PX4_ISFINITE(_pos_sp(i))) {
			valid = valid && PX4_ISFINITE(_pos(i));
		}

		if (PX4_ISFINITE(_vel_sp(i))) {
			valid = valid && PX4_ISFINITE(_vel(i)) && PX4_ISFINITE(_vel_dot(i));
		}
	}

	return valid;
}

void PositionControl::getLocalPositionSetpoint(vehicle_local_position_setpoint_s &local_position_setpoint) const
{
	local_position_setpoint.x = _pos_sp(0);
	local_position_setpoint.y = _pos_sp(1);
	local_position_setpoint.z = _pos_sp(2);
	local_position_setpoint.yaw = _yaw_sp;
	local_position_setpoint.yawspeed = _yawspeed_sp;
	local_position_setpoint.vx = _vel_sp(0);
	local_position_setpoint.vy = _vel_sp(1);
	local_position_setpoint.vz = _vel_sp(2);
	_acc_sp.copyTo(local_position_setpoint.acceleration);
	_thr_sp.copyTo(local_position_setpoint.thrust);
}

void PositionControl::getAttitudeSetpoint(vehicle_attitude_setpoint_s &attitude_setpoint) const
{
	ControlMath::thrustToAttitude(_thr_sp, _yaw_sp, attitude_setpoint);
	attitude_setpoint.yaw_sp_move_rate = _yawspeed_sp;
}
