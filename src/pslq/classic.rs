//! Textbook Ferguson–Bailey PSLQ, everything at full MPFR precision.

use super::state::State;
use super::{Outcome, PslqParams, RelationFinder};
use rug::Float;

pub struct ClassicPslq;

impl RelationFinder for ClassicPslq {
    fn name(&self) -> &'static str {
        "classic"
    }

    fn find(&self, x: &[Float], p: &PslqParams) -> Outcome {
        let mut st = State::new(x, p);
        let mut iterations = 0u64;
        loop {
            if let Some(outcome) = st.check(iterations) {
                return outcome;
            }
            iterations += 1;
            st.step();
        }
    }
}
