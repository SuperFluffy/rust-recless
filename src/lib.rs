extern crate ndarray;

use ndarray::prelude::*;
use ndarray::Data;
use ndarray::linalg::{
    general_mat_mul,
    general_mat_vec_mul,
};

/// The parameters of recursive least squares algorithm.
///
/// This struct contains all parameters involved in a recursive
/// least squares algorithm with exponential forgetting. For a review see
/// [Haykin's Adaptive Filter Theory][http://www.isbnsearch.org/isbn/9780132671453].
/// The implementation here does not implicitly take time into account. By making a choice of the
/// forgetting factor λ < 1 and shifting down old values of the input vector manually, the user can
/// get this algorithm to behave accordingly.
pub struct Rls<F> {

    /// The inverse forgetting factor λ^{-1}.
    inv_forgetting_factor: F,

    /// The gain vector used during the update of the inverse correlation matrix P(i) and the
    /// (tap) weight vector w(i).
    gain: Array1<F>,

    /// The *inverse correlation matrix*, initialized as P = δ^{-1} · 𝟙, with 𝟙 the unit matrix and
    /// δ an initialization factor recommended to be less than ~ 0.01 σ^2, the variance of the data
    /// sample to be analyzed (cf. [Haykin's Adaptive Filter Theory]).
    inverse_correlation: Array2<F>,

    /// The *(tap) weight vector* w(i) at time i used to produce the transversal filter's output
    /// y(i) = w(i) · u(i), where u(i) is the *(tap) input vector* at time i.
    weight: Array1<F>,

    /// The prior error used during the udpate of the inverse correlation matrix P(i) and the
    /// (tap) weight vector w(i). The prior error is calculated as the difference between the
    /// desired output and the filter output before an update.
    prior_error: F,

    // Two scratch matrices used for for update of the inverse correlation matrix.
    temp_a: Array2<F>,
    temp_b: Array2<F>,
}

impl<F: NdFloat> Rls<F> {
    /// Constructs a new Rls object with initialization factor δ and a weight vector of length n.
    pub fn new(initialization_factor: F, forgetting_factor: F, n: usize) -> Self {
        let weight = Array1::zeros(n);

        Rls::with_weight(initialization_factor, forgetting_factor, weight)
    }

    /// Constructs a new Rls object with initialization factor δ and pre-defined weight w.
    pub fn with_weight(initialization_factor: F, forgetting_factor: F, weight: Array1<F>) -> Self {
        let one = F::one();
        let zero = F::zero();

        let n = weight.len();

        let inv_forgetting_factor = one / forgetting_factor;

        let gain = Array1::zeros(n);
        let prior_error = zero;

        let inverse_correlation = Array2::from_elem([n,n], one/initialization_factor);

        let temp_a = Array2::zeros([n,n]);
        let temp_b = Array2::zeros([n,n]);

        Rls {
            inv_forgetting_factor,
            gain,
            inverse_correlation,
            weight,
            prior_error,
            temp_a,
            temp_b,
        }
    }

    /// Performs a recursive update of inverse correlation matrix and weight vector.
    pub fn update<S>(&mut self, input: &ArrayBase<S, Ix1>, target: F)
        where S: Data<Elem = F>
    {
        let one = F::one();
        let zero = F::zero();

        // Update the gain vector.
        general_mat_vec_mul(
            one,
            &self.inverse_correlation,
            input,
            one,
            &mut self.gain
        );
        let c = self.inv_forgetting_factor + self.gain.dot(&self.gain);
        self.gain /= c;

        // Calculate the prior error using the not yet updated tap weight.
        self.prior_error = target - self.weight.dot(&input);

        // Update the tap weight.
        self.weight.scaled_add(self.prior_error, &self.gain);

        // Update the inverse correlation matrix.
        for (mut row,a) in self.temp_a.genrows_mut().into_iter().zip(self.gain.iter()) {
            row.assign(input);
            row *= *a;
        }

        general_mat_mul(one, &self.temp_a, &self.inverse_correlation, zero, &mut self.temp_b);
        self.inverse_correlation -= &self.temp_b;
        self.inverse_correlation *= self.inv_forgetting_factor;
    }
}

impl<T> Rls<T> {
    pub fn gain_ref(&self) -> &Array1<T> {
        &self.gain
    }

    pub fn inverse_correlation_ref(&self) -> &Array2<T> {
        &self.inverse_correlation
    }

    pub fn inv_forgetting_factor_ref(&self) -> &T {
        &self.inv_forgetting_factor
    }

    pub fn weight_ref(&self) -> &Array1<T> {
        &self.weight
    }

    pub fn prior_error_ref(&self) -> &T {
        &self.prior_error
    }
}
