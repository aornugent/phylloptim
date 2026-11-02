// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_COST_CURVES_HPP_
#define PHYLLOPTIM_LEAF_COST_CURVES_HPP_

// The cost curves: each one's definition, `Curve<K>`, then its cost, profit and
// lambda, and the per-curve tables the solves dispatch on. A new curve is an
// enumerator, an arm in with_curve() in leaf.hpp, and a `Curve<K>` here, which
// `CurveDefinition` checks is complete.
//
// Part of `Leaf`: include <phylloptim/leaf_model.hpp>, which includes every
// part in order.

#include <phylloptim/leaf.hpp>

#include <string_view>

namespace phylloptim {

// --- a curve is one definition -----------------------------------------------
//
// `Curve<K>` states everything about curve K that the tables read: its name, the
// parameters it owns, its properties and its expressions. Every property is
// stated in every specialisation, with no base and no default: a defaulted
// `cost_reads_upstream = false` would solve the wrong first-order condition with
// no error. The concept checks each is present; review holds that none comes
// from a base.

// A parameter a curve owns: its slot in the pack, its bound field, the member it
// is, and its default, NaN where the curve refuses to solve until it is set.
struct CurveParameter {
  int slot;
  const char* field;
  double Leaf::*member;
  double default_value;
};

// What closed_form::optimise_into does with a curve: solves it; refuses it as
// unwritten (the expansion exists, its Newton slope does not); or refuses it as
// impossible (under the Log link lambda carries A, a fixed point).
enum class ClosedForm { Implemented, Unwritten, Impossible };

template <class T> inline constexpr bool is_curve_parameter_array = false;
template <std::size_t N>
inline constexpr bool
    is_curve_parameter_array<const std::array<CurveParameter, N>> = true;

// The properties, each read by the table of the same name or the branch that
// needs it, then the expressions at double, which the solves call:
//
//   link                          the benefit link composed with A
//   cost_reads_upstream           the cost reads the upstream potential too
//   cost_domain_ends_at_psi_crit  its own form ends there, not at the curve's end
//   objective_is_product          A times a gain factor that vanishes at psi_crit
//   cost_reads_leaf_temperature   the cost moves with Tleaf, so the energy
//                                 balance adds a cost-side slope
//   prepares_normalisers          normalised by quantities seeded before a search
//                                 (ProfitMax's |A|max and k_span), and with
//                                 unusable ones it has no objective
//   closed_form                   see ClosedForm
//   cost                          a difference curve's cost, in carbon, written
//                                 into hydraulic_cost_ (a gain or normalised
//                                 curve's is in its profit)
//   cost_slope                    dC/dpsi_stem, read by cost_deriv
//   cost_slope_upstream           dC/dpsi_upstream, where the cost reads it
//   lambda                        dC/dE at a stem potential, read by lambda_for
//   lambda_collar                 dC/dE on the collar route at the solved point
//   less_cost_temperature_slope   the EB term less the cost's slope in Tleaf,
//                                 where the cost reads it
//   profit                        the bound method that is its objective
template <class C>
concept CurveDefinition = requires(Leaf& l, const Leaf& cl, double x) {
  requires std::same_as<decltype(C::name), const std::string_view>;
  requires is_curve_parameter_array<decltype(C::parameters)>;
  requires std::same_as<decltype(C::link), const Leaf::BenefitLink>;
  requires std::same_as<decltype(C::cost_reads_upstream), const bool>;
  requires std::same_as<decltype(C::cost_domain_ends_at_psi_crit), const bool>;
  requires std::same_as<decltype(C::objective_is_product), const bool>;
  requires std::same_as<decltype(C::cost_reads_leaf_temperature), const bool>;
  requires std::same_as<decltype(C::prepares_normalisers), const bool>;
  requires std::same_as<decltype(C::closed_form), const ClosedForm>;
  { C::check_parameters(cl) } -> std::same_as<void>;
  { C::cost_slope(l, x, x) } -> std::same_as<double>;
  { C::lambda(l, x, x) } -> std::same_as<double>;
  { C::lambda_collar(l) } -> std::same_as<double>;
  requires std::same_as<decltype(C::profit),
                        double (Leaf::*const)(double, double)>;
} && (C::link != Leaf::BenefitLink::Identity || requires(Leaf& l, double x) {
  { C::cost(l, x, x) } -> std::same_as<double>;
  // Both ways for the slopes: one its property does not ask for is never called.
}) && (C::cost_reads_upstream == requires(Leaf& l, double x) {
  { C::cost_slope_upstream(l, x, x) } -> std::same_as<double>;
}) && (C::cost_reads_leaf_temperature ==
       requires(Leaf& l, double x, const Leaf::TemperaturePath& p) {
  { C::less_cost_temperature_slope(l, x, p) } -> std::same_as<double>;
});

// Aligned defaults: each level parameter (`CF77_lambda_`, `JS22_gamma`, `CMax_a`
// with `CMax_b` 0, `LeastCost_beta`) defaults to where its stem optimum has
// TF24's conductance at one reference: single-potential supply, psi_soil 1 MPa,
// PPFD 900, VPD 2 kPa, ca 40 Pa, leaf at 25 C, energy balance off, other traits
// at default. R's `align_curves()` solves them and `test-align-curves.R` fails
// when a change moves an optimum there: re-align every copy rather than loosen it.

template <> struct Curve<Leaf::CostCurve::TF24> {
  static constexpr std::string_view name = "TF24";
  // Traits, which TF24_floor's hydraulic half reads too.
  static constexpr std::array parameters{
      CurveParameter{par_TF24_beta2, "TF24_beta2", &Leaf::TF24_beta2, 1.5},
      CurveParameter{par_TF24_cost_scale, "TF24_cost_scale",
                     &Leaf::TF24_cost_scale, 7.5}};
  static constexpr Leaf::BenefitLink link = Leaf::BenefitLink::Identity;
  static constexpr bool cost_reads_upstream = false;
  static constexpr bool cost_domain_ends_at_psi_crit = false;
  static constexpr bool objective_is_product = false;
  static constexpr bool cost_reads_leaf_temperature = false;
  static constexpr bool prepares_normalisers = false;
  static constexpr ClosedForm closed_form = ClosedForm::Implemented;
  static void check_parameters(const Leaf&) {}

  static double cost(Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_upstream;
    l.hydraulic_cost_ = l.hydraulic_cost_TF_kernel(psi_stem);
    return l.hydraulic_cost_;
  }
  // The only cost slope that is a tangent of its own cost kernel.
  static double cost_slope(Leaf& l, double psi_stem, double psi_upstream) {
    using AD = xad::fwd<double>::active_type;
    (void)psi_upstream;
    AD ps_ad = psi_stem;
    xad::derivative(ps_ad) = 1.0;
    return xad::derivative(l.hydraulic_cost_TF_kernel(ps_ad));
  }
  // f cancels analytically, so this does not divide by the conductivity.
  static double lambda(const Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_upstream;
    const double f = l.proportion_of_conductivity(psi_stem);
    return l.TF24_cost_scale * l.TF24_beta2 * (l.stem_c / l.stem_b) * pow(psi_stem / l.stem_b, l.stem_c - 1.0) *
           pow(1.0 - f, l.TF24_beta2 - 1.0) / l.leaf_specific_conductance_max_;
  }
  static double lambda_collar(Leaf& l) {
    return l.marginal_cost_water_multilayer();
  }
  static constexpr auto profit = &Leaf::profit_psi_stem_TF;
};

template <> struct Curve<Leaf::CostCurve::CF77> {
  static constexpr std::string_view name = "CF77";
  static constexpr std::array parameters{CurveParameter{
      par_CF77_lambda, "CF77_lambda_", &Leaf::CF77_lambda_, 107750.8846}};
  static constexpr Leaf::BenefitLink link = Leaf::BenefitLink::Identity;
  static constexpr bool cost_reads_upstream = true;
  // A constant price need not overtake A'(E) before the stem's ceiling.
  static constexpr bool cost_domain_ends_at_psi_crit = true;
  static constexpr bool objective_is_product = false;
  static constexpr bool cost_reads_leaf_temperature = false;
  static constexpr bool prepares_normalisers = false;
  static constexpr ClosedForm closed_form = ClosedForm::Implemented;
  static void check_parameters(const Leaf& l) {
    // An input with an aligned default, so non-finite means a caller cleared it:
    // a NaN objective's argmax is a property of the bracket and looks plausible.
    if (!std::isfinite(l.CF77_lambda_)) {
      util::stop("the CF77 cost curve needs CF77_lambda_ set: it is the "
                 "PRESCRIBED marginal value of water in umol C (kg H2O)^-1, and "
                 "it was set to " + util::to_string(l.CF77_lambda_) + ". Its "
                 "default is aligned to TF24 (see align_curves()); "
                 "set_physiology and set_traits never touch it.");
    }
    // Only with the option on, so a caller who never sets it is not refused for
    // soil values the constant never reads. beta <= 0 prices water at infinity
    // or below zero, which pays the leaf to transpire to the wet bound.
    if (l.CF77_soil_beta_) {
      if (!(l.theta_fc_ > l.theta_w_)) {
        util::stop("CF77_soil_beta_ needs theta_fc_ > theta_w_; got theta_fc_ = " +
                   util::to_string(l.theta_fc_) + ", theta_w_ = " +
                   util::to_string(l.theta_w_));
      }
      const double beta = (l.theta_ - l.theta_w_) / (l.theta_fc_ - l.theta_w_);
      if (!std::isfinite(beta) || beta <= 0.0) {
        util::stop("CF77_soil_beta_ is on and the soil-moisture factor beta = "
                   "(theta_ - theta_w_)/(theta_fc_ - theta_w_) is " +
                   util::to_string(beta) + ", which prices water at infinity or "
                   "below zero. theta_ must exceed theta_w_; got theta_ = " +
                   util::to_string(l.theta_) + ", theta_w_ = " +
                   util::to_string(l.theta_w_));
      }
    }
  }

  // `CF77_lambda_ * E`, a carbon flux (umol C per kg water times kg m^-2 s^-1).
  // Straight after set_leaf_states_rates_from_psi_stem the memo returns that
  // point's E.
  static double cost(Leaf& l, double psi_stem, double psi_upstream) {
    l.hydraulic_cost_ = l.cf77_price() * l.operating_point_flux(psi_stem, psi_upstream);
    return l.hydraulic_cost_;
  }
  // lambda * E, and dE/dpsi_stem is kmax times the conductivity fraction,
  // which stem_curve_integral_deriv returns.
  static double cost_slope(Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_upstream;
    return l.cf77_price() * l.leaf_specific_conductance_max_ *
           l.stem_curve_integral_deriv(psi_stem);
  }
  // lambda * kmax * (G(psi_stem) - G(psi_up)), so -lambda * kmax * G'(psi_up),
  // with G' the table's own slope: the solve ran on the table.
  static double cost_slope_upstream(Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_stem;
    return -(l.cf77_price() * l.leaf_specific_conductance_max_ *
             l.stem_curve_integral_deriv(psi_upstream));
  }
  // Exact by construction: the cost is lambda*E, so dC/dpsi over dE/dpsi is the
  // price, independent of both potentials (the soil factor keys on theta).
  static double lambda(const Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_stem;
    (void)psi_upstream;
    return l.cf77_price();
  }
  static double lambda_collar(Leaf& l) {
    return l.cf77_price();
  }
  static constexpr auto profit = &Leaf::profit_psi_stem_CF77;
};

template <> struct Curve<Leaf::CostCurve::JS22> {
  static constexpr std::string_view name = "JS22";
  static constexpr std::array parameters{CurveParameter{
      par_JS22_gamma, "JS22_gamma", &Leaf::JS22_gamma, 0.4717887793}};
  static constexpr Leaf::BenefitLink link = Leaf::BenefitLink::Identity;
  static constexpr bool cost_reads_upstream = true;
  static constexpr bool cost_domain_ends_at_psi_crit = false;
  static constexpr bool objective_is_product = false;
  static constexpr bool cost_reads_leaf_temperature = false;
  static constexpr bool prepares_normalisers = false;
  static constexpr ClosedForm closed_form = ClosedForm::Unwritten;
  static void check_parameters(const Leaf& l) {
    // A trait with a default, still checked: a NaN maximises a NaN objective.
    if (!std::isfinite(l.JS22_gamma) || l.JS22_gamma < 0.0) {
      util::stop("JS22_gamma must be a finite, non-negative hydraulic unit cost "
                 "in umol C m^-2 s^-1 MPa^-2; got " + util::to_string(l.JS22_gamma));
    }
  }

  // gamma * (dpsi)^2, in carbon units. `dpsi*dpsi`, not `pow`: the optimiser's
  // inner loop (hazard 5).
  static double cost(Leaf& l, double psi_stem, double psi_upstream) {
    const double dpsi = psi_stem - psi_upstream;
    l.hydraulic_cost_ = l.JS22_gamma * dpsi * dpsi;
    return l.hydraulic_cost_;
  }
  static double cost_slope(Leaf& l, double psi_stem, double psi_upstream) {
    return 2.0 * l.JS22_gamma * (psi_stem - psi_upstream);
  }
  static double cost_slope_upstream(Leaf& l, double psi_stem, double psi_upstream) {
    return -2.0 * l.JS22_gamma * (psi_stem - psi_upstream);
  }
  // dC/dpsi over dE/dpsi, which is kmax*f(psi). Unlike TF24's, f does not cancel,
  // so this divides by it: safe inside the feasible bracket, where f never
  // reaches zero.
  static double lambda(const Leaf& l, double psi_stem, double psi_upstream) {
    const double f = l.proportion_of_conductivity(psi_stem);
    return 2.0 * l.JS22_gamma * (psi_stem - psi_upstream) /
           (l.leaf_specific_conductance_max_ * f);
  }
  // The cost reads the collar itself, so its partial over S joins the series.
  static double lambda_collar(Leaf& l) {
    double S, series;
    if (!l.collar_series(S, series)) {
      return util::na_value;
    }
    return lambda(l, l.opt_psi_stem_, l.opt_root_psi_) * series +
           cost_slope_upstream(l, l.opt_psi_stem_, l.opt_root_psi_) / S;
  }
  static constexpr auto profit = &Leaf::profit_psi_stem_JS22;
};

template <> struct Curve<Leaf::CostCurve::CMax> {
  static constexpr std::string_view name = "CMax";
  static constexpr std::array parameters{
      CurveParameter{par_CMax_a, "CMax_a", &Leaf::CMax_a, 0.6385596982},
      CurveParameter{par_CMax_b, "CMax_b", &Leaf::CMax_b, 0.0}};
  static constexpr Leaf::BenefitLink link = Leaf::BenefitLink::Identity;
  static constexpr bool cost_reads_upstream = true;
  static constexpr bool cost_domain_ends_at_psi_crit = false;
  static constexpr bool objective_is_product = false;
  static constexpr bool cost_reads_leaf_temperature = false;
  static constexpr bool prepares_normalisers = false;
  static constexpr ClosedForm closed_form = ClosedForm::Unwritten;
  static void check_parameters(const Leaf& l) {
    // Only CMax_a gets the sign check: CMax_b is signed, negative in the paper's
    // own convention (see the member).
    if (!std::isfinite(l.CMax_a) || l.CMax_a < 0.0) {
      util::stop("CMax_a must be a finite, non-negative slope for the CMax "
                 "marginal cost in umol C m^-2 s^-1 MPa^-2; got " +
                 util::to_string(l.CMax_a));
    }
    if (!std::isfinite(l.CMax_b)) {
      util::stop("CMax_b must be finite; got " + util::to_string(l.CMax_b));
    }
  }

  // The integral of `CMax_a*psi + CMax_b` from psi_upstream, zero at no flow:
  //
  //     C = (psi - psi_up) * [ CMax_a*(psi + psi_up)/2 + CMax_b ]
  //
  // Factored: the difference of squares cancels as the potentials converge (off by
  // 1.2e-04 at a 1e-12 drop), which a gradient's step reaches near no flow. For the
  // same reason JS22, this family's (2 gamma, -2 gamma psi_up) member, keeps its own.
  static double cost(Leaf& l, double psi_stem, double psi_upstream) {
    l.hydraulic_cost_ = (psi_stem - psi_upstream) *
                        (0.5 * l.CMax_a * (psi_stem + psi_upstream) + l.CMax_b);
    return l.hydraulic_cost_;
  }
  // Linear in the absolute potential: the form Anderegg et al. state.
  static double cost_slope(Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_upstream;
    return l.CMax_a * psi_stem + l.CMax_b;
  }
  // d/du of (psi - u)*(a*(psi + u)/2 + b) is -(a*u + b): the integrand at the
  // lower limit, as for any integral of the price.
  static double cost_slope_upstream(Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_stem;
    return -(l.CMax_a * psi_upstream + l.CMax_b);
  }
  // The marginal cost does not vanish at the wet end unless CMax_b is zero, so,
  // as CF77's and TF24_floor's, it can be boundary-soil with a non-zero cost there.
  static double lambda(const Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_upstream;
    const double f = l.proportion_of_conductivity(psi_stem);
    return (l.CMax_a * psi_stem + l.CMax_b) /
           (l.leaf_specific_conductance_max_ * f);
  }
  // As JS22's: the cost reads the collar as its lower limit.
  static double lambda_collar(Leaf& l) {
    double S, series;
    if (!l.collar_series(S, series)) {
      return util::na_value;
    }
    return lambda(l, l.opt_psi_stem_, l.opt_root_psi_) * series +
           cost_slope_upstream(l, l.opt_psi_stem_, l.opt_root_psi_) / S;
  }
  static constexpr auto profit = &Leaf::profit_psi_stem_CMax;
};

template <> struct Curve<Leaf::CostCurve::SOX> {
  static constexpr std::string_view name = "SOX";
  static constexpr std::array<CurveParameter, 0> parameters{};
  static constexpr Leaf::BenefitLink link = Leaf::BenefitLink::Log;
  static constexpr bool cost_reads_upstream = false;
  // The gain factor reaches zero at psi_crit and is negative past it.
  static constexpr bool cost_domain_ends_at_psi_crit = true;
  static constexpr bool objective_is_product = true;
  static constexpr bool cost_reads_leaf_temperature = false;
  static constexpr bool prepares_normalisers = false;
  static constexpr ClosedForm closed_form = ClosedForm::Impossible;
  static void check_parameters(const Leaf&) {}

  // Under the log link the cost is `-log g`, so dC/dpsi is `-g'/g`; `g` is zero
  // only at psi_crit, which no optimiser returns.
  static double cost_slope(Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_upstream;
    return -l.sox_reduction_deriv(psi_stem) / l.sox_reduction(psi_stem);
  }
  // A product: at its optimum A'g + Ag' = 0, so lambda = A'/E' = -A*(g'/g)/(kmax*f).
  // It carries a factor of A, since a product prices water as a share of
  // assimilation (the-models.Rmd), and diverges at psi_crit, where g = 0: never
  // pinned there, and not a number to evaluate at the dry end.
  static double lambda(const Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_upstream;
    const double g = l.sox_reduction(psi_stem);
    const double gp = l.sox_reduction_deriv(psi_stem);
    const double f = l.proportion_of_conductivity(psi_stem);
    return -l.assim_colimited_ * (gp / g) / (l.leaf_specific_conductance_max_ * f);
  }
  static double lambda_collar(Leaf& l) {
    double S, series;
    if (!l.collar_series(S, series)) {
      return util::na_value;
    }
    return lambda(l, l.opt_psi_stem_, l.opt_root_psi_) * series;
  }
  static constexpr auto profit = &Leaf::profit_psi_stem_SOX;
};

template <> struct Curve<Leaf::CostCurve::JW26> {
  static constexpr std::string_view name = "JW26";
  static constexpr std::array<CurveParameter, 0> parameters{};
  static constexpr Leaf::BenefitLink link = Leaf::BenefitLink::Log;
  static constexpr bool cost_reads_upstream = false;
  // As SOX's.
  static constexpr bool cost_domain_ends_at_psi_crit = true;
  static constexpr bool objective_is_product = true;
  static constexpr bool cost_reads_leaf_temperature = false;
  static constexpr bool prepares_normalisers = false;
  static constexpr ClosedForm closed_form = ClosedForm::Impossible;
  static void check_parameters(const Leaf&) {}

  // SOX's `-g'/g`, which collapses: with `g = 1 - psi/psi_crit` it is exactly
  // `1/(psi_crit - psi)`.
  static double cost_slope(Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_upstream;
    return 1.0 / (l.psi_crit - psi_stem);
  }
  // lambda = A / [(psi_crit - psi) kmax f(psi)], as the-models.Rmd quotes; it
  // carries A and diverges at psi_crit, as SOX's does.
  static double lambda(const Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_upstream;
    const double f = l.proportion_of_conductivity(psi_stem);
    return l.assim_colimited_ /
           ((l.psi_crit - psi_stem) * l.leaf_specific_conductance_max_ * f);
  }
  static double lambda_collar(Leaf& l) {
    double S, series;
    if (!l.collar_series(S, series)) {
      return util::na_value;
    }
    return lambda(l, l.opt_psi_stem_, l.opt_root_psi_) * series;
  }
  static constexpr auto profit = &Leaf::profit_psi_stem_JW26;
};

template <> struct Curve<Leaf::CostCurve::ProfitMax> {
  static constexpr std::string_view name = "ProfitMax";
  static constexpr std::array<CurveParameter, 0> parameters{};
  static constexpr Leaf::BenefitLink link = Leaf::BenefitLink::Scaled;
  static constexpr bool cost_reads_upstream = false;
  static constexpr bool cost_domain_ends_at_psi_crit = false;
  static constexpr bool objective_is_product = false;
  // Its thermal cost, with the energy balance on.
  static constexpr bool cost_reads_leaf_temperature = true;
  static constexpr bool prepares_normalisers = true;
  static constexpr ClosedForm closed_form = ClosedForm::Unwritten;
  static void check_parameters(const Leaf&) {}

  // C = [k(psi_s) - k(psi)]/k_span, so dC/dpsi = kmax*|f'|/k_span. The thermal
  // half's slope is less_cost_temperature_slope.
  static double cost_slope(Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_upstream;
    const double f = l.proportion_of_conductivity(psi_stem);
    const double abs_fprime =
        f * (l.stem_c / l.stem_b) * pow(psi_stem / l.stem_b, l.stem_c - 1.0);
    return l.leaf_specific_conductance_max_ * abs_fprime / l.profitmax_k_span_;
  }
  // Cost (k_soil - kmax f)/k_span + TC(Tleaf), in carbon units by |A|max:
  //
  //     lambda = |A|max * [ |f'|/(f*k_span) + (dTC/dT)(dT/dE) ]
  //
  // The thermal term vanishes with either gate off.
  static double lambda(Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_upstream;
    const double f = l.proportion_of_conductivity(psi_stem);
    if (!(f > 0.0) || !(l.profitmax_k_span_ > 0.0) || !(l.profitmax_A_max_ > 0.0)) {
      return util::na_value;
    }
    // |f'| = (c/b)(psi/b)^(c-1) f, positive because f is decreasing.
    const double abs_fprime =
        (l.stem_c / l.stem_b) * pow(psi_stem / l.stem_b, l.stem_c - 1.0) * f;

    double thermal = 0.0;
    if (l.use_thermal_cost_ && l.use_energy_balance_) {
      const double E = l.operating_point_flux(psi_stem, l.supply_psi_soil_scalar());
      double dT_dE = 0.0;
      const double Tleaf = l.leaf_temp_from_E(E, &dT_dE);
      // The same slope less_cost_temperature_slope subtracts from dJ/dpsi, so the
      // reported lambda and the first-order condition cannot disagree about it.
      thermal = l.thermal_cost_temp_deriv(Tleaf) * dT_dE;
    }
    return l.profitmax_A_max_ *
           (abs_fprime / (f * l.profitmax_k_span_) + thermal);
  }
  // The conductance term reads psi_stem alone, so it takes the series factor;
  // the thermal term, where the gates put one in, reads Tleaf and so E alone.
  static double lambda_collar(Leaf& l) {
    double S, series;
    if (!l.collar_series(S, series)) {
      return util::na_value;
    }
    const double f = l.proportion_of_conductivity(l.opt_psi_stem_);
    if (!(f > 0.0) || !(l.profitmax_k_span_ > 0.0) || !(l.profitmax_A_max_ > 0.0)) {
      return util::na_value;
    }
    const double abs_fprime = (l.stem_c / l.stem_b) *
                              pow(l.opt_psi_stem_ / l.stem_b, l.stem_c - 1.0) * f;
    double thermal = 0.0;
    if (l.use_thermal_cost_ && l.use_energy_balance_) {
      double dT_dE = 0.0;
      const double Tleaf = l.leaf_temp_from_E(
          l.operating_point_flux(l.opt_psi_stem_, l.opt_root_psi_), &dT_dE);
      thermal = l.thermal_cost_temp_deriv(Tleaf) * dT_dE;
    }
    return l.profitmax_A_max_ *
           (abs_fprime / (f * l.profitmax_k_span_) * series + thermal);
  }
  // dTC/dpsi = TC'(Tleaf) * tau, the slope lambda reports. One expression with
  // the subtraction: returned alone, the product is rounded before it is
  // subtracted, and the energy-balance rows move.
  static double less_cost_temperature_slope(const Leaf& l, double term,
                                            const Leaf::TemperaturePath& p) {
    return term - l.thermal_cost_temp_deriv(p.Tleaf) * p.tau;
  }
  static constexpr auto profit = &Leaf::profit_psi_stem_ProfitMax;
};

template <> struct Curve<Leaf::CostCurve::TF24_floor> {
  static constexpr std::string_view name = "TF24_floor";
  // Its own price; its hydraulic half reads TF24's parameters.
  static constexpr std::array parameters{
      CurveParameter{par_TF24_floor_lambda_o, "TF24_floor_lambda_o",
                     &Leaf::TF24_floor_lambda_o, util::na_value}};
  static constexpr Leaf::BenefitLink link = Leaf::BenefitLink::Identity;
  static constexpr bool cost_reads_upstream = true;
  static constexpr bool cost_domain_ends_at_psi_crit = false;
  static constexpr bool objective_is_product = false;
  static constexpr bool cost_reads_leaf_temperature = false;
  static constexpr bool prepares_normalisers = false;
  static constexpr ClosedForm closed_form = ClosedForm::Implemented;
  static void check_parameters(const Leaf& l) {
    // No default, unlike CF77_lambda_: at lambda_o = 0 this curve IS TF24, so an
    // unset price would run TF24 under this curve's name. An explicit zero is
    // accepted (the reduction test's case); the hydraulic half is TF24's traits.
    if (!std::isfinite(l.TF24_floor_lambda_o) || l.TF24_floor_lambda_o < 0.0) {
      util::stop("the TF24_floor cost curve needs TF24_floor_lambda_o set: it is the "
                 "PRESCRIBED marginal cost of water as transpiration goes to "
                 "zero, in umol C (kg H2O)^-1, non-negative, and NA until you "
                 "assign one -- set_physiology and set_traits never touch it. It "
                 "has no default because 0 would silently make this curve TF24. "
                 "Set `$TF24_floor_lambda_o` on the leaf, or pass "
                 "`TF24_floor_lambda_o =` to leaf_solve() / leaf_batch(); got " +
                 util::to_string(l.TF24_floor_lambda_o));
    }
  }

  // TF24's cost plus CF77's price, both carbon fluxes. Each term here and in the
  // slopes and lambda is its parent's own expression in its operand order, so
  // lambda_o = 0 and TF24_cost_scale = 0 reduce bit-exactly, under FMA contraction
  // too (hazards.md, "One solver: the long form"). Reassociating breaks that;
  // test_tf24_floor_reduces guards it.
  static double cost(Leaf& l, double psi_stem, double psi_upstream) {
    l.hydraulic_cost_ = l.hydraulic_cost_TF_kernel(psi_stem) +
                        l.TF24_floor_lambda_o * l.operating_point_flux(psi_stem, psi_upstream);
    return l.hydraulic_cost_;
  }
  static double cost_slope(Leaf& l, double psi_stem, double psi_upstream) {
    using AD = xad::fwd<double>::active_type;
    (void)psi_upstream;
    AD ps_ad = psi_stem;
    xad::derivative(ps_ad) = 1.0;
    return xad::derivative(l.hydraulic_cost_TF_kernel(ps_ad)) +
           l.TF24_floor_lambda_o * l.leaf_specific_conductance_max_ *
               l.stem_curve_integral_deriv(psi_stem);
  }
  // TF24's half reads psi_stem alone; the floor's lambda_o*E is CF77's arm.
  static double cost_slope_upstream(Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_stem;
    return -(l.TF24_floor_lambda_o * l.leaf_specific_conductance_max_ *
             l.stem_curve_integral_deriv(psi_upstream));
  }
  // TF24's lambda plus lambda_o. Called, not repeated: at lambda_o == 0 it is
  // TF24's bit for bit, since that ends in a division no FMA can absorb.
  static double lambda(const Leaf& l, double psi_stem, double psi_upstream) {
    (void)psi_upstream;
    return l.lambda_TF24(psi_stem) + l.TF24_floor_lambda_o;
  }
  // TF24's part takes the series factor; the price reads E alone.
  static double lambda_collar(Leaf& l) {
    double S, series;
    if (!l.collar_series(S, series)) {
      return util::na_value;
    }
    return l.lambda_TF24(l.opt_psi_stem_) * series + l.TF24_floor_lambda_o;
  }
  static constexpr auto profit = &Leaf::profit_psi_stem_TF24_floor;
};

template <> struct Curve<Leaf::CostCurve::LeastCost> {
  static constexpr std::string_view name = "LeastCost";
  static constexpr std::array parameters{CurveParameter{
      par_LeastCost_beta, "LeastCost_beta", &Leaf::LeastCost_beta, 51.59667042}};
  static constexpr Leaf::BenefitLink link = Leaf::BenefitLink::Log;
  // `log W` reads upstream through E, the only log-link cost that does.
  static constexpr bool cost_reads_upstream = true;
  // It sits at psi_crit on most rows, dark and light.
  static constexpr bool cost_domain_ends_at_psi_crit = true;
  static constexpr bool objective_is_product = false;
  // W's capacity term reads vcmax_ at Tleaf.
  static constexpr bool cost_reads_leaf_temperature = true;
  static constexpr bool prepares_normalisers = false;
  static constexpr ClosedForm closed_form = ClosedForm::Impossible;
  static void check_parameters(const Leaf& l) {
    // An aligned default, refused when cleared for CF77's reason. Strictly
    // positive, unlike TF24_floor_lambda_o: zero leaves `A/E`, which rises without
    // bound as E -> 0 and puts every optimum at the wet bound.
    if (!std::isfinite(l.LeastCost_beta) || l.LeastCost_beta <= 0.0) {
      util::stop("the LeastCost cost curve needs LeastCost_beta set: it is "
                 "Prentice et al. (2014)'s beta, the unit cost of maintaining "
                 "carboxylation capacity relative to that of transpiration, "
                 "DIMENSIONLESS (both fluxes are weighed in mol m^-2 s^-1), "
                 "strictly positive; its default is aligned to TF24 (see "
                 "align_curves()), and set_physiology and set_traits never touch "
                 "it. The published value is 146. Set `$LeastCost_beta` on the "
                 "leaf, or pass "
                 "`LeastCost_beta =` to leaf_solve() / leaf_batch(); got " +
                 util::to_string(l.LeastCost_beta));
    }
  }

  // Under the log link the cost is `+log W`, `W = E_mol + beta*Vcmax_mol`
  // (`leastcost_flux`), so dC/dpsi is `W'/W` with `W' = kg_to_mol_h2o * kmax * f`:
  // the opposite sign to SOX's `-g'/g`, as W rises where g falls. Vcmax's slope
  // in Tleaf is less_cost_temperature_slope.
  static double cost_slope(Leaf& l, double psi_stem, double psi_upstream) {
    return kg_to_mol_h2o * l.leaf_specific_conductance_max_ *
           l.stem_curve_integral_deriv(psi_stem) /
           l.leastcost_flux(psi_stem, psi_upstream);
  }
  // `(dW/du)/W` with `dW/du = -kg_to_mol_h2o * kmax * G'(u)`: CF77's arm with
  // kg-to-mol in place of the price, over cost_slope's `W`.
  static double cost_slope_upstream(Leaf& l, double psi_stem, double psi_upstream) {
    return -(kg_to_mol_h2o * l.leaf_specific_conductance_max_ *
             l.stem_curve_integral_deriv(psi_upstream)) /
           l.leastcost_flux(psi_stem, psi_upstream);
  }
  // `lambda = A*kg_to_mol_h2o/W`, reading the solved `assim_colimited_` as SOX's
  // does: valid only at the potential last evaluated.
  static double lambda(Leaf& l, double psi_stem, double psi_upstream) {
    if (!l.use_energy_balance_) {
      return l.assim_colimited_ * kg_to_mol_h2o /
             l.leastcost_flux(psi_stem, psi_upstream);
    }
    // W's capacity term moves with E through Tleaf, adding beta*Vcmax_mol'(T)*dT/dE:
    // the slope less_cost_temperature_slope takes in the first-order condition.
    double dT_dE = 0.0;
    const double Tleaf =
        l.leaf_temp_from_E(l.operating_point_flux(psi_stem, psi_upstream), &dT_dE);
    double dW_dE = kg_to_mol_h2o;
    if (dT_dE != 0.0) {
      dW_dE += l.LeastCost_beta * umol_to_mol * l.vcmax_temp_deriv(Tleaf) * dT_dE;
    }
    return l.assim_colimited_ * dW_dE / l.leastcost_flux(psi_stem, psi_upstream);
  }
  // The cost reads E alone, so its lambda carries over unchanged.
  static double lambda_collar(Leaf& l) {
    return lambda(l, l.opt_psi_stem_, l.opt_root_psi_);
  }
  // C = log W with `vcmax_` at Tleaf, so d(log W)/dpsi gains
  // beta*umol_to_mol*dVcmax/dT*tau/W (cost_slope holds only the E half). W reads
  // the total flux, its stomatal part from gc.
  static double less_cost_temperature_slope(const Leaf& l, double term,
                                            const Leaf::TemperaturePath& p) {
    const double dVcmax_dT = (p.vcmax_up - p.vcmax_dn) / (2.0 * p.h);
    const double W = p.gc_s / p.gc_const * kg_to_mol_h2o +
                     l.LeastCost_beta * l.vcmax_ * umol_to_mol +
                     p.r_star * kg_to_mol_h2o;
    return term - l.LeastCost_beta * umol_to_mol * dVcmax_dT * p.tau / W;
  }
  static constexpr auto profit = &Leaf::profit_psi_stem_LeastCost;
};

// Every enumerator has a complete definition. Asserted per curve, so a missing
// piece is reported against the curve and the requirement it fails.
template <Leaf::CostCurve K>
constexpr bool curve_is_defined() {
  static_assert(CurveDefinition<Curve<K>>,
                "this curve's Curve<K> fails CurveDefinition; the notes name the requirement");
  // The normalisers the solves seed are ProfitMax's (prepare_profitmax), so a
  // second normalised curve needs its own before it can say true.
  static_assert(!Curve<K>::prepares_normalisers || K == Leaf::CostCurve::ProfitMax,
                "prepares_normalisers seeds ProfitMax's normalisers");
  // ClosedForm says why a closed form is impossible: the Log link.
  static_assert((Curve<K>::closed_form == ClosedForm::Impossible) ==
                    (Curve<K>::link == Leaf::BenefitLink::Log),
                "a closed form is impossible exactly for a Log-link curve");
  return true;
}

template <std::size_t... I>
constexpr bool every_curve_is_defined(std::index_sequence<I...>) {
  return (curve_is_defined<static_cast<Leaf::CostCurve>(I)>() && ...);
}

static_assert(every_curve_is_defined(
    std::make_index_sequence<static_cast<std::size_t>(Leaf::n_cost_curves)>{}));

// Every curve's name, in enum order.
template <std::size_t... I>
constexpr std::array<std::string_view, sizeof...(I)> curve_names_of(
    std::index_sequence<I...>) {
  return {Curve<static_cast<Leaf::CostCurve>(I)>::name...};
}

inline constexpr auto curve_names = curve_names_of(
    std::make_index_sequence<static_cast<std::size_t>(Leaf::n_cost_curves)>{});

// set_model() finds a curve by its name and R a parameter by its slot, so
// neither may repeat.
constexpr bool curve_names_are_unique() {
  for (std::size_t i = 0; i < curve_names.size(); ++i) {
    for (std::size_t j = i + 1; j < curve_names.size(); ++j) {
      if (curve_names[i] == curve_names[j]) return false;
    }
  }
  return true;
}
static_assert(curve_names_are_unique(), "two curves share a name");

template <std::size_t... I>
constexpr bool curve_slots_are_unique(std::index_sequence<I...>) {
  std::array<bool, n_leaf_pars> seen{};
  bool unique = true;
  const auto mark = [&](const auto& parameters) {
    for (const CurveParameter& p : parameters) {
      if (seen[static_cast<std::size_t>(p.slot)]) unique = false;
      seen[static_cast<std::size_t>(p.slot)] = true;
    }
  };
  (mark(Curve<static_cast<Leaf::CostCurve>(I)>::parameters), ...);
  return unique;
}
static_assert(curve_slots_are_unique(std::make_index_sequence<
                  static_cast<std::size_t>(Leaf::n_cost_curves)>{}),
              "two curve parameters share a slot");

template <Leaf::CostCurve K, int Slot>
constexpr double Leaf::curve_default() {
  constexpr const auto& ps = Curve<K>::parameters;
  constexpr auto it = std::find_if(
      ps.begin(), ps.end(), [](const CurveParameter& p) { return p.slot == Slot; });
  static_assert(it != ps.end(), "curve_default: the curve owns no parameter in this slot");
  return it->default_value;
}

// dC/dpsi_stem for one cost curve: the only place a cost curve enters a
// derivative, shared by both routes (the collar route multiplies it by
// dpsi_stem/dpsi and adds cost_deriv_upstream<K>).
template <Leaf::CostCurve K>
inline double Leaf::cost_deriv(double psi_stem, double psi_upstream) {
  return Curve<K>::cost_slope(*this, psi_stem, psi_upstream);
}

// Whether a curve's cost reads the upstream potential as well as psi_stem, which
// matters on the collar route, where upstream is the decision variable. A curve
// that reads it but says false solves the wrong first-order condition with no
// error; `test_collar_foc_matches_profit` checks every curve.
template <Leaf::CostCurve K>
constexpr bool Leaf::cost_reads_upstream() {
  return Curve<K>::cost_reads_upstream;
}

// Whether a curve's own form ends at psi_crit. Where it does not, the marginal
// cost diverges as the stem dries, so the dry end is the end of the stem curve
// (vulnerability_psi_end).
template <Leaf::CostCurve K>
constexpr bool Leaf::cost_domain_ends_at_psi_crit() {
  return Curve<K>::cost_domain_ends_at_psi_crit;
}

// SOX and JW26 maximise A g, g falling to 0 at psi_crit. With A <= 0 at every
// opening the supremum is 0 as g -> 0, driving a dark leaf to psi_crit, so a
// leaf with A <= 0 at the dry end (A rises with opening) is shut at the wet
// bound instead. The argmax jumps at the light compensation point regardless.
template <Leaf::CostCurve K>
constexpr bool Leaf::objective_is_product() {
  return Curve<K>::objective_is_product;
}

template <Leaf::CostCurve K>
inline bool Leaf::dark_under_product(double psi_dry, double psi_upstream) {
  static_assert(objective_is_product<K>(),
                "dark_under_product is only defined for a product objective");
  set_leaf_states_rates_from_psi_stem(psi_dry, psi_upstream);
  return !(assim_colimited_ > 0.0);
}

template <Leaf::CostCurve K>
inline double Leaf::dry_end() const {
  if constexpr (cost_domain_ends_at_psi_crit<K>()) {
    return psi_crit;
  } else {
    return vulnerability_psi_end(stem_b, stem_c);
  }
}

// dC/dpsi_upstream at fixed psi_stem, for the curves `cost_reads_upstream` names.
// The collar FOC subtracts it with a unit chain factor, since upstream is the
// collar there.
template <Leaf::CostCurve K>
inline double Leaf::cost_deriv_upstream(double psi_stem, double psi_upstream) {
  static_assert(cost_reads_upstream<K>(),
                "cost_deriv_upstream is only defined where the cost reads upstream");
  return Curve<K>::cost_slope_upstream(*this, psi_stem, psi_upstream);
}

// --- Marginal cost of water -------------------------------------------------
//
// The bound per-curve lambdas, each its curve's definition's.

inline double Leaf::lambda_TF24(double psi_stem) const {
  return Curve<CostCurve::TF24>::lambda(*this, psi_stem, util::na_value);
}

inline double Leaf::lambda_JS22(double psi_stem, double psi_upstream) const {
  return Curve<CostCurve::JS22>::lambda(*this, psi_stem, psi_upstream);
}

inline double Leaf::lambda_CMax(double psi_stem, double psi_upstream) const {
  return Curve<CostCurve::CMax>::lambda(*this, psi_stem, psi_upstream);
}

inline double Leaf::lambda_TF24_floor(double psi_stem) const {
  return Curve<CostCurve::TF24_floor>::lambda(*this, psi_stem, util::na_value);
}

inline double Leaf::lambda_SOX(double psi_stem, double psi_upstream) const {
  return Curve<CostCurve::SOX>::lambda(*this, psi_stem, psi_upstream);
}

inline double Leaf::lambda_JW26(double psi_stem, double psi_upstream) const {
  return Curve<CostCurve::JW26>::lambda(*this, psi_stem, psi_upstream);
}

inline double Leaf::lambda_ProfitMax(double psi_stem) {
  return Curve<CostCurve::ProfitMax>::lambda(*this, psi_stem, util::na_value);
}

inline double Leaf::marginal_cost_water() const {
  return lambda_TF24(opt_psi_stem_);
}

inline double Leaf::marginal_cost_water_molar() const {
  // umol CO2 (kg H2O)^-1 -> mol CO2 (mol H2O)^-1
  return marginal_cost_water() * umol_to_mol / kg_to_mol_h2o;
}

inline double Leaf::marginal_cost_water_multilayer() {
  // S is the soil->collar conductance, positive by construction (hazard 2).
  const double S = dE_from_soil_dpsi_collar(opt_root_psi_, supply_psi_soil());
  if (!std::isfinite(S) || S <= 0.0) {
    return util::na_value;
  }
  // f is the STEM vulnerability curve at the collar suction: kmax*f(psi_r) is
  // dE_stem/dpsi_r, the stem-side conductance where the two paths meet.
  const double f_r = proportion_of_conductivity(opt_root_psi_);
  return lambda_TF24(opt_psi_stem_) *
         (1.0 + leaf_specific_conductance_max_ * f_r / S);
}

// Pure: no write to hydraulic_cost_, so the AD pass cannot scribble model state
// while probing. The caching is the double entry point's job.
template <typename T, typename U>
inline T Leaf::hydraulic_cost_TF_kernel(const T& psi_stem, const U& scale,
                                        const U& b, const U& c,
                                        const U& beta2) const {
  return scale *
         pow((1 - proportion_of_conductivity_kernel<T, U>(psi_stem, b, c)), beta2);
}

template <typename T>
inline T Leaf::hydraulic_cost_TF_kernel(T psi_stem) const {
  return hydraulic_cost_TF_kernel<T, double>(psi_stem, TF24_cost_scale, stem_b,
                                             stem_c, TF24_beta2);
}

template <typename T>
inline T Leaf::cost_slope_at(const T& sigma, const leaf_pars<T>& pars) const {
  // Named rather than passed inline: a row borrows its input by reference, and a
  // temporary is not a thing to hand one.
  const T b =
      phylloptim::weibull_b_from_P50<T>(pars[par_stem_P50], pars[par_stem_c]);
  return kernel_slope_at<T>(
      [this]<class U>(const U& s, const U& scale, const U& bb, const U& c,
                      const U& beta2) -> U {
        return hydraulic_cost_TF_kernel<U, U>(s, scale, bb, c, beta2);
      },
      "hydraulic cost in stem potential", sigma, pars[par_TF24_cost_scale], b,
      pars[par_stem_c], pars[par_TF24_beta2]);
}

// The curve itself is closed form, so the (P50, c) chain comes out of the one
// derivation rather than a supplied row: only the pre-integrated G(psi) is tabulated and
// needs its rows supplied (stem_integral_at).
template <typename T>
inline T Leaf::hydraulic_cost_TF_kernel(const T& psi_stem,
                                        const leaf_pars<T>& pars) const {
  return hydraulic_cost_TF_kernel<T>(
      psi_stem, pars[par_TF24_cost_scale],
      phylloptim::weibull_b_from_P50<T>(pars[par_stem_P50], pars[par_stem_c]),
      pars[par_stem_c], pars[par_TF24_beta2]);
}

// The bound difference-curve costs, each its curve's definition's.
inline double Leaf::hydraulic_cost_TF(double psi_stem) {
  return Curve<CostCurve::TF24>::cost(*this, psi_stem, util::na_value);
}

// Sicangco et al. (2026) Eqns 8-9: the fraction of maximum PSII damage at this
// leaf temperature, in [0,1]; `r` follows gsthermal, not the printed Eqn 8 (see
// use_thermal_cost_). Exactly 0.0 with the gate off, so callers add it freely.
inline double Leaf::thermal_cost_at(double leaf_temp) const {
  if (!use_thermal_cost_) {
    return 0.0;
  }
  const double span = T50_ - Tcrit_;
  if (!(span > 0.0)) {
    util::stop("thermal cost needs T50 > Tcrit; got T50 = " +
               util::format_double(T50_) + ", Tcrit = " +
               util::format_double(Tcrit_));
  }
  return 1.0 / (1.0 + std::exp(-(2.0 / span) * (leaf_temp - T50_)));
}

inline double Leaf::thermal_cost_temp_deriv(double leaf_temp) const {
  if (!use_thermal_cost_) {
    return 0.0;
  }
  const double span = T50_ - Tcrit_;
  const double TC = thermal_cost_at(leaf_temp);
  return (2.0 / span) * TC * (1.0 - TC);
}

// Profit functions

inline double Leaf::profit_psi_stem_TF(double psi_stem, double psi_upstream) {
set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_TF(psi_stem);

  return benefit_ - cost;
}

inline double Leaf::hydraulic_cost_CF77(double psi_stem, double psi_upstream) {
  return Curve<CostCurve::CF77>::cost(*this, psi_stem, psi_upstream);
}

// The CF77 price with the optional Medlyn soil-moisture factor, rising without
// bound toward the wilting point:
//
//   beta = (theta - theta_w) / (theta_fc - theta_w),   price = lambda / beta
//
// beta is clamped at 1, so CF77_lambda_ is a floor and theta >= theta_fc is plain
// CF77 exactly. The off path returns the member, not a division by 1.0, so it
// stays bit-identical to plain CF77.
inline double Leaf::cf77_price() const {
  if (!CF77_soil_beta_) { return CF77_lambda_; }
  const double beta = (theta_ - theta_w_) / (theta_fc_ - theta_w_);
  return CF77_lambda_ / (beta > 1.0 ? 1.0 : beta);
}

inline double Leaf::profit_psi_stem_CF77(double psi_stem,
                                                 double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_CF77(psi_stem, psi_upstream);

  return benefit_ - cost;
}

inline double Leaf::hydraulic_cost_JS22(double psi_stem, double psi_upstream) {
  return Curve<CostCurve::JS22>::cost(*this, psi_stem, psi_upstream);
}

inline double Leaf::profit_psi_stem_JS22(double psi_stem,
                                         double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_JS22(psi_stem, psi_upstream);

  return benefit_ - cost;
}

inline double Leaf::hydraulic_cost_CMax(double psi_stem, double psi_upstream) {
  return Curve<CostCurve::CMax>::cost(*this, psi_stem, psi_upstream);
}

inline double Leaf::profit_psi_stem_CMax(double psi_stem,
                                         double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_CMax(psi_stem, psi_upstream);

  return benefit_ - cost;
}

inline double Leaf::hydraulic_cost_TF24_floor(double psi_stem, double psi_upstream) {
  return Curve<CostCurve::TF24_floor>::cost(*this, psi_stem, psi_upstream);
}

// See the declaration for why zero on the other curves is a statement about the
// curve and not about the cost.
inline double Leaf::shadow_cost() const {
  return cost_curve_ == CostCurve::TF24_floor
           ? TF24_floor_lambda_o * transpiration_
           : 0.0;
}

inline double Leaf::profit_psi_stem_TF24_floor(double psi_stem,
                                            double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_TF24_floor(psi_stem, psi_upstream);

  return benefit_ - cost;
}

// --- the product-objective family --------------------------------------------
//
// Eller's reduction factor: the conductivity fraction rescaled to 1 at zero
// tension and exactly 0 at `psi_crit`, through the k_crit_fraction psi_crit and
// Sperry's k_crit also read, so the three cannot disagree.
// Clamped at zero because the no-flow branch evaluates past psi_crit, where g < 0
// and A = -R_d made a shut leaf's product positive. `g < 0`, not `std::max`, so a
// NaN propagates.
inline double Leaf::sox_reduction(double psi_stem) const {
  const double g = (proportion_of_conductivity(psi_stem) - k_crit_fraction) /
                   (1.0 - k_crit_fraction);
  return g < 0.0 ? 0.0 : g;
}

// Its slope, negative: `|f'| = f*(c/b)*(psi/b)^(c-1)`, lambda_TF24's analytic form.
inline double Leaf::sox_reduction_deriv(double psi_stem) const {
  const double f = proportion_of_conductivity(psi_stem);
  const double abs_fprime =
      f * (stem_c / stem_b) * pow(psi_stem / stem_b, stem_c - 1.0);
  return -abs_fprime / (1.0 - k_crit_fraction);
}

// A product, not in carbon units: the argmax is log A + log g's, but the value
// must not be compared with a TF24 profit.
// ⚠️ hydraulic_cost_ is poisoned (hazard 8): with no subtracted cost, leaving it
// would report whichever curve ran last on a reused Leaf.
inline double Leaf::profit_psi_stem_SOX(double psi_stem,
                                        double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  hydraulic_cost_ = util::na_value;
  return assim_colimited_ * sox_reduction(psi_stem);
}

// Jones' linear reduction factor on the derived `psi_crit`, clamped at zero, NaN
// propagating, as `sox_reduction` is.
inline double Leaf::jw26_reduction(double psi_stem) const {
  const double g = 1.0 - psi_stem / psi_crit;
  return g < 0.0 ? 0.0 : g;
}

// Constant, the whole difference from SOX, whose slope follows the vulnerability
// curve. Zero once the clamp binds.
inline double Leaf::jw26_reduction_deriv(double psi_stem) const {
  return psi_stem > psi_crit ? 0.0 : -1.0 / psi_crit;
}

inline double Leaf::profit_psi_stem_JW26(double psi_stem,
                                         double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  hydraulic_cost_ = util::na_value;
  return assim_colimited_ * jw26_reduction(psi_stem);
}

// --- least-cost / cost-per-unit-carbon ---------------------------------------

// `W = E_mol + beta*Vcmax_mol`, mol m^-2 s^-1; see the declaration for why a
// dimensionless beta needs both conversions. It reads `vcmax_`, not `vcmax_25`, so
// with the energy balance on it moves with the solved leaf temperature, as the
// ratio Prentice et al. minimise; cost_deriv leaves that motion out and
// dprofit_energy_balance_term adds it.
inline double Leaf::leastcost_flux(double psi_stem, double psi_upstream) {
  return operating_point_flux(psi_stem, psi_upstream) * kg_to_mol_h2o +
         LeastCost_beta * vcmax_ * umol_to_mol;
}

// A ratio, not a carbon flux; `hydraulic_cost_` is poisoned as for SOX (hazard 8).
inline double Leaf::profit_psi_stem_LeastCost(double psi_stem,
                                              double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  hydraulic_cost_ = util::na_value;
  return assim_colimited_ / leastcost_flux(psi_stem, psi_upstream);
}

inline double Leaf::lambda_LeastCost(double psi_stem, double psi_upstream) {
  return Curve<CostCurve::LeastCost>::lambda(*this, psi_stem, psi_upstream);
}

// ===========================================================================
// Sperry et al. (2017) ProfitMax, as Sicangco et al. (2026) implement it
// ---------------------------------------------------------------------------
// Sperry maximises `Profit = CG - HC`, both normalised. Times |A|max that is
// `A - lambda*(k(psi_soil)-k(psi))` with the same argmax, but lambda* moves with
// every driver, so a prescribed lambda is not this model:
//
//     CG      = A(psi)/|A|max
//     HC      = [k(psi_soil)-k(psi)] / [k(psi_soil)-kcrit]
//     lambda* = |A|max / [k(psi_soil) - kcrit]
//
// The profit is in Sperry's units, so a figure overlays the paper's. HC runs 0 to
// 1 however wide the vulnerability curve, where the TF24 cost becomes negligible
// as stem_b widens. kcrit = k_crit_fraction * kmax, so psi_crit is P95, as Sperry.
// HC measures loss from psi_soil on one soil potential, and on a layered profile
// from the collar at which uptake is zero (root_zero_E, prepare_collar_solve's
// wet bracket end).
inline double Leaf::profitmax_reference_psi() {
  if (supply_is_single_layer()) {
    return supply_psi_soil_scalar();
  }
  double bound_a = 0.0, bound_b = 0.0;
  if (!prepare_collar_solve<CostCurve::TF24>(bound_a, bound_b)) {
    return util::na_value;
  }
  return bound_a;
}

inline void Leaf::prepare_profitmax_norms() {
  // A non-finite reference is propagated, not refused: this runs on leaves whose
  // drivers are not seated yet (the batch route seeds |A|max first). The
  // optimiser's `!(profitmax_k_span_ > 0)` guard is the one place that decides it.
  const double psi_soil = profitmax_reference_psi();
  profitmax_k_soil_ =
      leaf_specific_conductance_max_ * proportion_of_conductivity(psi_soil);
  profitmax_k_span_ =
      profitmax_k_soil_ -
      leaf_specific_conductance_max_ * proportion_of_conductivity(psi_crit);
}

inline void Leaf::prepare_profitmax_at(double A_max) {
  prepare_profitmax_norms();
  if (!(std::isfinite(A_max) && A_max > 0.0)) {
    util::stop("prepare_profitmax_at needs a finite positive |A|max -- it is the "
               "normaliser held fixed across a perturbation, so an unset one "
               "means the base point was never solved; got " +
               util::to_string(A_max));
  }
  profitmax_A_max_ = A_max;
}

inline void Leaf::prepare_profitmax() {
  const double psi_soil = profitmax_reference_psi();
  prepare_profitmax_norms();

  // |A|max by the one solver, not a scan, whose argmax is piecewise constant in
  // the parameters (hazard 3). One maximisation of A suffices: |A| peaks at a
  // minimum of A only where A < 0, at closure or along a whole stream that cannot
  // fix carbon (below). dA/dpsi is dJ/dpsi plus the cost derivative, for an
  // identity link. Closure is excluded, as gsthermal drops it (`A[E == 0] <- NA`):
  // the search starts one probe inside the wet bound.
  const auto A_of = [&](double psi) {
    set_leaf_states_rates_from_psi_stem(psi, psi_soil);
    return assim_colimited_;
  };
  const auto dA_of = [&](double psi, bool* ok) {
    const double dJ = dprofit_dpsi_stem<CostCurve::TF24>(psi, psi_soil, ok);
    return dJ + cost_deriv<CostCurve::TF24>(psi, psi_soil);
  };
  const double psi_wet = stem_wet_bound(psi_soil);
  const double width = psi_crit - psi_wet;
  const double lo = psi_wet + 1e-9 * width;   // the open end: closure excluded
  double a_hi = -std::numeric_limits<double>::infinity();
  if (width > 0.0) {
    util::maximise_over_closed_interval_foc(A_of, dA_of, lo, psi_crit, 0,
                                            collar_root_tol,
                                            static_cast<size_t>(ci_niter), &a_hi);
  }
  double a_max = (std::isfinite(a_hi) && a_hi > 0.0) ? a_hi : 0.0;
  // The reachable case where max|A| is not max A: A negative along the whole
  // stream, a leaf too hot to cover respiration. Normalising by zero would leave
  // no objective, so maximise -A, and only then.
  if (!(a_max > 0.0) && width > 0.0) {
    double neg = -std::numeric_limits<double>::infinity();
    util::maximise_over_closed_interval_foc(
        [&](double psi) { return -A_of(psi); },
        [&](double psi, bool* ok) { return -dA_of(psi, ok); }, lo, psi_crit, 0,
        collar_root_tol, static_cast<size_t>(ci_niter), &neg);
    if (std::isfinite(neg) && neg > 0.0) a_max = neg;
  }
  profitmax_A_max_ = a_max;
}

// CG - (HC + TC) at one candidate potential, with the normalisers prepare_profitmax
// seeded. Writes carbon_gain_, hydraulic_cost_norm_ and thermal_cost_ so a caller
// plotting the paper's Figure 2 can read the three components off the object.
inline double Leaf::profit_psi_stem_ProfitMax(double psi_stem,
                                              double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  carbon_gain_ = (profitmax_A_max_ > 0.0)
                     ? assim_colimited_ / profitmax_A_max_
                     : 0.0;
  hydraulic_cost_norm_ =
      (profitmax_k_soil_ -
       leaf_specific_conductance_max_ * proportion_of_conductivity(psi_stem)) /
      profitmax_k_span_;
  // At the operating point's leaf temperature; with the energy balance off TC is
  // constant along the stream, as in Sicangco's ProfitMaxTC without it.
  thermal_cost_ = thermal_cost_at(use_energy_balance_
                                      ? leaf_temp_from_E(transpiration_)
                                      : leaf_temp_);

  return carbon_gain_ - (hydraulic_cost_norm_ + thermal_cost_);
}

inline std::vector<double> Leaf::profitmax_curve(int n) {
  if (n < 2) {
    util::stop("profitmax_curve needs at least 2 points");
  }
  prepare_profitmax();
  const double psi_soil = supply_psi_soil_scalar();
  const double psi_wet = stem_wet_bound(psi_soil);
  const double step = (psi_crit - psi_wet) / double(n - 1);
  std::vector<double> out(static_cast<std::size_t>(5 * n));
  for (int i = 0; i < n; ++i) {
    const double p = psi_wet + step * double(i);
    const double profit = profit_psi_stem_ProfitMax(p, psi_soil);
    const std::size_t k = static_cast<std::size_t>(i);
    const std::size_t N = static_cast<std::size_t>(n);
    out[k] = p;
    out[N + k] = carbon_gain_;
    out[2 * N + k] = hydraulic_cost_norm_;
    out[3 * N + k] = thermal_cost_;
    out[4 * N + k] = profit;
  }
  return out;
}

inline void Leaf::set_profitmax_degenerate_state(double psi_stem,
                                                 double psi_upstream) {
  set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);
  profit_ = 0;
  carbon_gain_ = 0;
  hydraulic_cost_norm_ = 0;
  thermal_cost_ = 0;
  lambda_emergent_ = util::na_value;
}

// --- runtime curve selection --------------------------------------------------

inline std::string Leaf::curve_name(int curve) {
  if (curve < 0 || curve >= n_cost_curves) {
    return "unknown";
  }
  return std::string(curve_names[static_cast<std::size_t>(curve)]);
}

// Which link each curve composes with assimilation. SOX and JW26 are products
// `A*g`, LeastCost a ratio `A/W`: the link is the same `h = log`, and the sign
// (-log g against +log W) lands in `cost_deriv`.
template <Leaf::CostCurve K>
constexpr Leaf::BenefitLink Leaf::benefit_link() {
  return Curve<K>::link;
}

// `h'(A)`, all the derivative needs from the link. The Log arm differentiates
// `log A + log g`, which shares `A*g`'s argmax. `log A` needs `A > 0`, and at
// full closure `A = -R_d`: a non-positive `A` returns NaN, which the caller's
// feasibility flag carries outward.
template <Leaf::CostCurve K>
inline double Leaf::benefit_link_deriv(double A) const {
  if constexpr (benefit_link<K>() == BenefitLink::Log) {
    return A > 0.0 ? 1.0 / A : util::na_value;
  } else if constexpr (benefit_link<K>() == BenefitLink::Scaled) {
    return profitmax_A_max_ > 0.0 ? 1.0 / profitmax_A_max_
                                   : util::na_value;
  } else {
    (void)A;
    return 1.0;
  }
}

// What must be set before the search: each curve's own checks.
template <Leaf::CostCurve K>
inline void Leaf::check_cost_parameters() {
  Curve<K>::check_parameters(*this);
}

// --- the per-curve tables ---------------------------------------------------
//
// Each reads its curve's definition, which CurveDefinition holds complete.

// Which objective. Not all one kind of number: five are carbon, two products
// (`profit_psi_stem_SOX`) and one a ratio (`profit_psi_stem_LeastCost`).
template <Leaf::CostCurve K>
inline double Leaf::profit_psi_stem_for(double psi_stem, double psi_upstream) {
  return (this->*Curve<K>::profit)(psi_stem, psi_upstream);
}

// Which lambda to report, normalised to two arguments (TF24's reads psi_stem
// alone).
template <Leaf::CostCurve K>
inline double Leaf::lambda_for(double psi_stem, double psi_upstream) {
  return Curve<K>::lambda(*this, psi_stem, psi_upstream);
}

// lambda_emergent on the COLLAR route: dC/dE along the path the solve moves on,
// the collar psi_r free and psi_stem following. Uptake gives dE/dpsi_r = S and the
// stem path dpsi_stem/dpsi_r = (S + kmax f(psi_r)) / (kmax f(psi_stem)), so the
// part of a cost reading psi_stem alone takes the series correction
//
//     lambda = lambda_stem * (1 + kmax f(psi_r) / S),
//
// the part reading only E carries over unchanged, and a cost that also reads the
// collar (JS22's drop, CMax's lower limit) adds its partial over S. Each curve's
// definition writes its own; each must equal dA/dE at an interior optimum:
// test_collar_lambda_emergent_is_dA_dE.
template <Leaf::CostCurve K>
inline double Leaf::lambda_collar_for() {
  return Curve<K>::lambda_collar(*this);
}

inline bool Leaf::collar_series(double& S, double& series) {
  S = dE_from_soil_dpsi_collar(opt_root_psi_, supply_psi_soil());
  if (!std::isfinite(S) || S <= 0.0) {
    return false;
  }
  series = 1.0 + leaf_specific_conductance_max_ *
                     proportion_of_conductivity(opt_root_psi_) / S;
  return true;
}

} // namespace phylloptim

#endif
