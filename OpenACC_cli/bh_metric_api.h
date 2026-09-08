#ifndef BH_METRIC_API_H
#define BH_METRIC_API_H

/**
 * Contract shared by every generated metric header.
 *
 * A generated header `bh_<name>_metric.h` opens a namespace `bh_<name>` and
 * provides, all as `template <class FP>`:
 *
 *   metric(p, x, g)                g_{mu nu} at x
 *   geodesic_rhs(p, x, v, acc)     acc^mu = -Gamma^mu_{ab} v^a v^b
 *   observer(p, x, u)              the camera's 4-velocity, up to scale
 *   event_fields(p, x, f)          one scalar per termination event
 *   event_guards(p, x, gd)         guard scalars, concatenated
 *
 * plus the N_PARAM / N_EVENT / EVENT_KIND tables declared below.  `p` is a flat
 * parameter vector holding both the metric's own parameters and the scene's;
 * PARAM_NAME and PARAM_DEFAULT describe it to the CLI.
 *
 * Nothing here assumes a (t, r, theta, phi) chart.  A ray stops when an event
 * fires, and an event is a sign condition on a scalar field, which is a
 * statement that survives any change of coordinates.  A Cartesian chart
 * describes its sky sphere as `R_sky - sqrt(x*x + y*y + z*z)` and its disk as
 * the sign change of `z`; a Boyer-Lindquist chart writes the same two events as
 * `R_sky - r` and the sign change of `theta - pi/2`.  The tracer does not need
 * to know which it is looking at.
 */

//! How an event's scalar field turns into "the ray stops here".
enum bh_event_kind
{
    //! Fires while f(x) < 0: a region the ray falls into, e.g. inside a horizon.
    EVENT_REGION_NEG = 0,
    //! Fires when f(x) goes from positive to negative across a step.
    EVENT_CROSS_POS_TO_NEG,
    //! Fires when f(x) goes from negative to positive across a step.
    EVENT_CROSS_NEG_TO_POS
};

#endif // BH_METRIC_API_H
