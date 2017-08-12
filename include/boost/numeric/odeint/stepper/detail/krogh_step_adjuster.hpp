/*
 boost/numeric/odeint/stepper/detail/pid_step_adjuster.hpp

 [begin_description]
 Implementation of the stepsize controller for the controlled adams bashforth moulton stepper.
 [end_description]

 Copyright 2017 Valentin Noah Hartmann

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_KROGH_STEP_ADJUSTER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_KROGH_STEP_ADJUSTER_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>

#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>

#include <cmath>

#include <boost/array.hpp>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template<
class State,
class Value = double,
class Deriv = State,
class Time = double,
class Algebra = typename algebra_dispatcher< State >::algebra_type,
class Operations = typename operations_dispatcher< Deriv >::operations_type
>
struct krogh_step_adjuster
{
public:
    static double threshold() { return 0.9; };

    typedef State state_type;
    typedef Value value_type;
    typedef Deriv deriv_type;
    typedef Time time_type;

    typedef Algebra algebra_type;
    typedef Operations operations_type;

    typedef rotating_buffer<state_type, 3> error_storage_type;
    typedef rotating_buffer<time_type, 3> time_storage_type;

    krogh_step_adjuster(double abs_tol = 1e-6, double rel_tol = 1e-6, time_type dtmax = 1.0)
    :m_dtmax(dtmax), m_init(0),
    m_abs_tol(abs_tol), m_rel_tol(rel_tol)
    {};

    time_type adjust_stepsize(const size_t steps, time_type dt, state_type &err, const state_type &x, const deriv_type &dxdt)
    {
        using std::abs;
        m_algebra.for_each3( err , x , dxdt ,
                typename operations_type::template rel_error< value_type >( m_abs_tol , m_rel_tol , 1.0 , 1.0 * abs(get_unit_value( dt )) ) );

        m_phi[0] = err;

        // actual operation
        value_type w = 0.9;

        for(size_t j=0; j<err.size(); ++j)
        {
            m_phi[0][j] = log(err[j]) - steps*log(dt);
        }

        value_type ratio;

        if(m_init == 0)
        {
            for(size_t j=0; j<err.size(); ++j)
            {
                err[j] = pow(err[j], 1/(value_type(steps)));
            }
            ratio = 1/m_algebra.norm_inf(err);
        }
        else
        {
            if(m_init == 1)
            {
                // initialize
                m_r[0] = err;
                m_r[1] = err;
                m_r[2] = err;

                for(size_t j=0; j<err.size(); ++j)
                {
                    m_r[0][j] = (w*m_phi[0][j]+(1-2*w)*m_phi[1][j])/((1-w)*(1-w));
                    m_r[1][j] = (2*w*m_phi[0][j]+(1-3*w)*m_phi[1][j])/((1-w)*(1-w)*(1-w));
                    m_r[2][j] = (3*w*m_phi[0][j]+(1-4*w)*m_phi[1][j])/((1-w)*(1-w)*(1-w)*(1-w));
                }
            }

            else
            {
                for(size_t j=0; j<err.size(); ++j)
                {
                    m_r[0][j] = m_phi[0][j]+w*m_r[0][j];
                    m_r[1][j] = m_r[0][j]+w*m_r[1][j];
                    m_r[2][j] = m_r[1][j]+w*m_r[2][j];
                }
            }

            state_type a = err;
            for(size_t j=0; j<err.size(); ++j)
            {
                a[j] = (1-w*w)/w*m_r[0][j]-((1-w)*(1-w))/w*m_r[1][j];
                // a[j] = (1-w)/(w*w)*((1+w+w*w)*m_r[0][j] + (-2+w+w*w)*m_r[1][j]+(1-2*w+w*w)*m_r[2][j]);
                err[j] = exp(a[j]/steps);
            }

            ratio = 1/(m_algebra.norm_inf(err))/dt;
        }

        value_type kappa = 1.0;
        ratio = 1.0 + kappa*atan((ratio - 1) / kappa);

        if(ratio*dt >= m_dtmax)
        {
            ratio = m_dtmax / dt;
        }

        if(ratio >= threshold())
        {
            m_phi.rotate();

            ++m_init;
        }
        else
        {
            m_init = 0;
        }

        return dt * static_cast<time_type>(ratio);
    };

private:
    algebra_type m_algebra;

    time_type m_dtmax;

    boost::array<state_type, 3> m_r;

    size_t m_init;
    double m_abs_tol;
    double m_rel_tol;

    error_storage_type m_phi;
};

} // detail
} // odeint
} // numeric
} // boost

#endif