/// Bond types in a molecular graph.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BondType {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl BondType {
    /// Bond order × 10 (aromatic = 15 for integer math).
    pub fn order_x10(&self) -> u8 {
        match self {
            BondType::Single => 10,
            BondType::Double => 20,
            BondType::Triple => 30,
            BondType::Aromatic => 15,
        }
    }

    /// Contribution to total valence for implicit H calculation.
    pub fn valence_contribution(&self) -> u8 {
        match self {
            BondType::Single => 1,
            BondType::Double => 2,
            BondType::Triple => 3,
            BondType::Aromatic => 1,
        }
    }
}
