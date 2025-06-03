/****************************************************************************
 *
 *   Copyright (c) 2019-2023 PX4 Development Team. All rights reserved.
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
 * @file RateControl.cpp
 */

#include "rate_control.hpp"
#include <px4_platform_common/defines.h>
#include <cmath>

using namespace matrix;

float s_roll = 0.0f;
float s_pitch = 0.0f;
float s_yaw = 0.0f;

float tau_eq_roll = 0.0f;
float tau_eq_pitch = 0.0f;
float tau_eq_yaw = 0.0f;

float tau_sw_roll = 0.0f;
float tau_sw_pitch = 0.0f;
float tau_sw_yaw = 0.0f;

float tauR = 0.0f;
float tauP = 0.0f;
float tauY = 0.0f;

void RateControl::setPidGains(const Vector3f &P, const Vector3f &I, const Vector3f &D)
{
	_gain_p = P;
	_gain_i = I;
	_gain_d = D;
}

void RateControl::setSmcGains(const matrix::Vector3f &I, const matrix::Vector3f &K1, const matrix::Vector3f &K2, const matrix::Vector3f &c1, const matrix::Vector3f &c2, const matrix::Vector3f &eta)
{
	_gain_I = I;
	_gain_k1 = K1;
	_gain_k2 = K2;
	_gain_c1 = c1;
	_gain_c2 = c2;
	_gain_eta = eta;
}

void RateControl::setSaturationStatus(const Vector3<bool> &saturation_positive,
				      const Vector3<bool> &saturation_negative)
{
	_control_allocator_saturation_positive = saturation_positive;
	_control_allocator_saturation_negative = saturation_negative;
}

void RateControl::setPositiveSaturationFlag(size_t axis, bool is_saturated)
{
	if (axis < 3) {
		_control_allocator_saturation_positive(axis) = is_saturated;
	}
}

void RateControl::setNegativeSaturationFlag(size_t axis, bool is_saturated)
{
	if (axis < 3) {
		_control_allocator_saturation_negative(axis) = is_saturated;
	}
}

// Vector3f RateControl::update(const Vector3f &att_cur, const Vector3f &att_sp, const Vector3f &rate, const Vector3f &rate_sp, const Vector3f &angular_accel,
			     	// const float dt, const bool landed)
// {
// 	// angular rates error
// 	Vector3f rate_error = rate_sp - rate;

// 	// PID control with feed forward
// 	const Vector3f torque = _gain_p.emult(rate_error) + _rate_int - _gain_d.emult(angular_accel) + _gain_ff.emult(rate_sp);

// 	// update integral only if we are not landed
// 	if (!landed) {
// 		updateIntegral(rate_error, dt);
// 	}

// 	return torque;
// }

/* Modified Version of Update */

Vector3f RateControl::update(const Vector3f &att_cur, const Vector3f &att_sp, const Vector3f &rate, const Vector3f &rate_sp, const Vector3f &angular_accel,
			     const float dt, const bool landed)
{
	// angular rates error
	Vector3f rate_error = rate_sp - rate;
	Vector3f att_err = att_sp - att_cur;

	// sliding surface
	s_roll = (1.f + _gain_c2(0)) * rate_error(0) + _gain_c1(0) * att_err(0);
	s_pitch = (1.f + _gain_c2(1)) * rate_error(1) + _gain_c1(1) * att_err(1);
	s_yaw = (1.f + _gain_c2(2)) * rate_error(2) + _gain_c1(2) * att_err(2);

	// Equivalent control
	tau_eq_roll = _gain_I(0) * (0.f + (_gain_c1(0) / (1.f + _gain_c2(0))) * rate_error(0));
	tau_eq_pitch = _gain_I(1) * (0.f + (_gain_c1(1) / (1.f + _gain_c2(1))) * rate_error(1));
	tau_eq_yaw = _gain_I(2) * (0.f + (_gain_c1(2) / (1.f + _gain_c2(2))) * rate_error(2));

	tau_sw_roll = _gain_I(0) * (_gain_k2(0) / (1.f + _gain_c2(0))) * s_roll + _gain_I(0) * (_gain_k1(0) / (1.f + _gain_c2(0))) * tanhf(_gain_eta(0) * s_roll);
	tau_sw_pitch = _gain_I(1) * (_gain_k2(1) / (1.f + _gain_c2(1))) * s_pitch + _gain_I(1) * (_gain_k1(1) / (1.f + _gain_c2(1))) * tanhf(_gain_eta(1) * s_pitch);
	tau_sw_yaw = _gain_I(2) * (_gain_k2(2) / (1.f + _gain_c2(2))) * s_yaw + _gain_I(2) * (_gain_k1(2) / (1.f + _gain_c2(2))) * tanhf(_gain_eta(2) * s_yaw);

	tauR = tau_eq_roll + tau_sw_roll;
	tauP = tau_eq_pitch + tau_sw_pitch;
	tauY = tau_eq_yaw + tau_sw_yaw;

	// SMC control with feed forward
	const Vector3f torque = {tauR, tauP, tauY};

	// update integral only if we are not landed
	if (!landed) {
		updateIntegral(rate_error, dt);
	}

	return torque;
}

// Vector3f RateControl::update(const Vector3f &rate, const Vector3f &rate_sp,
// 							 const Vector3f &angular_accel, const float dt, const bool landed)
// {
// 	// Compute rate setpoint derivative
// 	Vector3f rate_sp_dot;
// 	if (_first_update){
// 		rate_sp_dot = Vector3f(0.f, 0.f, 0.f);
// 		_first_update = false;
// 	} else{
// 		rate_sp_dot = (rate_sp - _rate_sp_prev) / dt;
// 	}
// 	_rate_sp_prev = rate_sp;

// 	// Define error
// 	Vector3f e = rate_sp - rate;
// 	// Compute error derivative
//     Vector3f e_dot = rate_sp_dot - angular_accel;
//     // Compute sliding surface s = e_dot + c * e
//     Vector3f s = e_dot + _gain_c * e;

// 	// Compute torque based on desired s_dot = -k2 * s - k1 * tanh(eta * s)
//     Vector3f torque;
//     for (int i = 0; i < 3; i++) {
//         float tanh_term = std::tanh(_gain_eta(i) * s(i));
//         // Since torque = I * angular_accel, and we want s_dot = -k2*s - k1*tanh(eta*s),
//         // we approximate torque = K * (rate_sp_dot - (k2*s + k1*tanh(eta*s)))
//         torque(i) = _gain_k(i) * (rate_sp_dot(i) - (_gain_k2(i) * s(i) + _gain_k1(i) * tanh_term));
//     }

//     return torque;
// }

/* The modified Program Ends Here */

void RateControl::updateIntegral(Vector3f &rate_error, const float dt)
{
	for (int i = 0; i < 3; i++) {
		// prevent further positive control saturation
		if (_control_allocator_saturation_positive(i)) {
			rate_error(i) = math::min(rate_error(i), 0.f);
		}

		// prevent further negative control saturation
		if (_control_allocator_saturation_negative(i)) {
			rate_error(i) = math::max(rate_error(i), 0.f);
		}

		// I term factor: reduce the I gain with increasing rate error.
		// This counteracts a non-linear effect where the integral builds up quickly upon a large setpoint
		// change (noticeable in a bounce-back effect after a flip).
		// The formula leads to a gradual decrease w/o steps, while only affecting the cases where it should:
		// with the parameter set to 400 degrees, up to 100 deg rate error, i_factor is almost 1 (having no effect),
		// and up to 200 deg error leads to <25% reduction of I.
		float i_factor = rate_error(i) / math::radians(400.f);
		i_factor = math::max(0.0f, 1.f - i_factor * i_factor);

		// Perform the integration using a first order method
		float rate_i = _rate_int(i) + i_factor * _gain_i(i) * rate_error(i) * dt;

		// do not propagate the result if out of range or invalid
		if (PX4_ISFINITE(rate_i)) {
			_rate_int(i) = math::constrain(rate_i, -_lim_int(i), _lim_int(i));
		}
	}
}

void RateControl::getRateControlStatus(rate_ctrl_status_s &rate_ctrl_status)
{
	rate_ctrl_status.rollspeed_integ = _rate_int(0);
	rate_ctrl_status.pitchspeed_integ = _rate_int(1);
	rate_ctrl_status.yawspeed_integ = _rate_int(2);
}
