use std::fmt;

/// Chemical elements relevant for drug-like molecules.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Element {
    H,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F,
    Ne,
    Na,
    Mg,
    Al,
    Si,
    P,
    S,
    Cl,
    Ar,
    K,
    Ca,
    Ti,
    V,
    Cr,
    Mn,
    Fe,
    Co,
    Ni,
    Cu,
    Zn,
    Ga,
    Ge,
    As,
    Se,
    Br,
    Kr,
    Rb,
    Sr,
    Zr,
    Mo,
    Ru,
    Rh,
    Pd,
    Ag,
    Cd,
    In,
    Sn,
    Sb,
    Te,
    I,
    Xe,
    Cs,
    Ba,
    Pt,
    Au,
    Hg,
    Tl,
    Pb,
    Bi,
    Unknown,
}

impl Element {
    /// Parse from element symbol string (case-sensitive: "C", "Cl", "Br").
    pub fn from_symbol(s: &str) -> Option<Self> {
        match s {
            "H" => Some(Element::H),
            "He" => Some(Element::He),
            "Li" => Some(Element::Li),
            "Be" => Some(Element::Be),
            "B" => Some(Element::B),
            "C" => Some(Element::C),
            "N" => Some(Element::N),
            "O" => Some(Element::O),
            "F" => Some(Element::F),
            "Ne" => Some(Element::Ne),
            "Na" => Some(Element::Na),
            "Mg" => Some(Element::Mg),
            "Al" => Some(Element::Al),
            "Si" => Some(Element::Si),
            "P" => Some(Element::P),
            "S" => Some(Element::S),
            "Cl" => Some(Element::Cl),
            "Ar" => Some(Element::Ar),
            "K" => Some(Element::K),
            "Ca" => Some(Element::Ca),
            "Ti" => Some(Element::Ti),
            "V" => Some(Element::V),
            "Cr" => Some(Element::Cr),
            "Mn" => Some(Element::Mn),
            "Fe" => Some(Element::Fe),
            "Co" => Some(Element::Co),
            "Ni" => Some(Element::Ni),
            "Cu" => Some(Element::Cu),
            "Zn" => Some(Element::Zn),
            "Ga" => Some(Element::Ga),
            "Ge" => Some(Element::Ge),
            "As" => Some(Element::As),
            "Se" => Some(Element::Se),
            "Br" => Some(Element::Br),
            "Kr" => Some(Element::Kr),
            "Rb" => Some(Element::Rb),
            "Sr" => Some(Element::Sr),
            "Zr" => Some(Element::Zr),
            "Mo" => Some(Element::Mo),
            "Ru" => Some(Element::Ru),
            "Rh" => Some(Element::Rh),
            "Pd" => Some(Element::Pd),
            "Ag" => Some(Element::Ag),
            "Cd" => Some(Element::Cd),
            "In" => Some(Element::In),
            "Sn" => Some(Element::Sn),
            "Sb" => Some(Element::Sb),
            "Te" => Some(Element::Te),
            "I" => Some(Element::I),
            "Xe" => Some(Element::Xe),
            "Cs" => Some(Element::Cs),
            "Ba" => Some(Element::Ba),
            "Pt" => Some(Element::Pt),
            "Au" => Some(Element::Au),
            "Hg" => Some(Element::Hg),
            "Tl" => Some(Element::Tl),
            "Pb" => Some(Element::Pb),
            "Bi" => Some(Element::Bi),
            _ => None,
        }
    }

    /// Atomic number.
    pub fn atomic_number(&self) -> u8 {
        match self {
            Element::H => 1,
            Element::He => 2,
            Element::Li => 3,
            Element::Be => 4,
            Element::B => 5,
            Element::C => 6,
            Element::N => 7,
            Element::O => 8,
            Element::F => 9,
            Element::Ne => 10,
            Element::Na => 11,
            Element::Mg => 12,
            Element::Al => 13,
            Element::Si => 14,
            Element::P => 15,
            Element::S => 16,
            Element::Cl => 17,
            Element::Ar => 18,
            Element::K => 19,
            Element::Ca => 20,
            Element::Ti => 22,
            Element::V => 23,
            Element::Cr => 24,
            Element::Mn => 25,
            Element::Fe => 26,
            Element::Co => 27,
            Element::Ni => 28,
            Element::Cu => 29,
            Element::Zn => 30,
            Element::Ga => 31,
            Element::Ge => 32,
            Element::As => 33,
            Element::Se => 34,
            Element::Br => 35,
            Element::Kr => 36,
            Element::Rb => 37,
            Element::Sr => 38,
            Element::Zr => 40,
            Element::Mo => 42,
            Element::Ru => 44,
            Element::Rh => 45,
            Element::Pd => 46,
            Element::Ag => 47,
            Element::Cd => 48,
            Element::In => 49,
            Element::Sn => 50,
            Element::Sb => 51,
            Element::Te => 52,
            Element::I => 53,
            Element::Xe => 54,
            Element::Cs => 55,
            Element::Ba => 56,
            Element::Pt => 78,
            Element::Au => 79,
            Element::Hg => 80,
            Element::Tl => 81,
            Element::Pb => 82,
            Element::Bi => 83,
            Element::Unknown => 0,
        }
    }

    /// Symbol string.
    pub fn symbol(&self) -> &'static str {
        match self {
            Element::H => "H",
            Element::He => "He",
            Element::Li => "Li",
            Element::Be => "Be",
            Element::B => "B",
            Element::C => "C",
            Element::N => "N",
            Element::O => "O",
            Element::F => "F",
            Element::Ne => "Ne",
            Element::Na => "Na",
            Element::Mg => "Mg",
            Element::Al => "Al",
            Element::Si => "Si",
            Element::P => "P",
            Element::S => "S",
            Element::Cl => "Cl",
            Element::Ar => "Ar",
            Element::K => "K",
            Element::Ca => "Ca",
            Element::Ti => "Ti",
            Element::V => "V",
            Element::Cr => "Cr",
            Element::Mn => "Mn",
            Element::Fe => "Fe",
            Element::Co => "Co",
            Element::Ni => "Ni",
            Element::Cu => "Cu",
            Element::Zn => "Zn",
            Element::Ga => "Ga",
            Element::Ge => "Ge",
            Element::As => "As",
            Element::Se => "Se",
            Element::Br => "Br",
            Element::Kr => "Kr",
            Element::Rb => "Rb",
            Element::Sr => "Sr",
            Element::Zr => "Zr",
            Element::Mo => "Mo",
            Element::Ru => "Ru",
            Element::Rh => "Rh",
            Element::Pd => "Pd",
            Element::Ag => "Ag",
            Element::Cd => "Cd",
            Element::In => "In",
            Element::Sn => "Sn",
            Element::Sb => "Sb",
            Element::Te => "Te",
            Element::I => "I",
            Element::Xe => "Xe",
            Element::Cs => "Cs",
            Element::Ba => "Ba",
            Element::Pt => "Pt",
            Element::Au => "Au",
            Element::Hg => "Hg",
            Element::Tl => "Tl",
            Element::Pb => "Pb",
            Element::Bi => "Bi",
            Element::Unknown => "*",
        }
    }

    /// Default valence(s) for implicit hydrogen calculation.
    pub fn default_valences(&self) -> &'static [u8] {
        match self {
            Element::B => &[3],
            Element::C => &[4],
            Element::N => &[3, 5],
            Element::O => &[2],
            Element::P => &[3, 5],
            Element::S => &[2, 4, 6],
            Element::F => &[1],
            Element::Cl => &[1, 3, 5, 7],
            Element::Br => &[1, 3, 5, 7],
            Element::I => &[1, 3, 5, 7],
            Element::Si => &[4],
            Element::Al => &[3],
            Element::Sn => &[2, 4],
            Element::Pb => &[2, 4],
            _ => &[0],
        }
    }

    /// Is this an organic subset atom?
    pub fn is_organic_subset(&self) -> bool {
        matches!(
            self,
            Element::B
                | Element::C
                | Element::N
                | Element::O
                | Element::P
                | Element::S
                | Element::F
                | Element::Cl
                | Element::Br
                | Element::I
        )
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

#[derive(Clone, Debug)]
pub struct Atom {
    pub element: Element,
    pub charge: i8,
    pub isotope: Option<u16>,
    /// Total H count (explicit + implicit); assigned after graph construction.
    pub h_count: u8,
    /// Explicit H count from bracket atom syntax, e.g. `[CH3]` → Some(3).
    /// None means "no bracket H specification" (compute implicitly).
    pub explicit_h: Option<u8>,
    pub aromatic: bool,
    pub in_ring: bool,
    pub ring_sizes: smallvec::SmallVec<[u8; 4]>,
    pub map_num: Option<u16>,
}

impl Atom {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            charge: 0,
            isotope: None,
            h_count: 0,
            explicit_h: None,
            aromatic: false,
            in_ring: false,
            ring_sizes: smallvec::SmallVec::new(),
            map_num: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_symbol() {
        assert_eq!(Element::from_symbol("C"), Some(Element::C));
        assert_eq!(Element::from_symbol("Cl"), Some(Element::Cl));
        assert_eq!(Element::from_symbol("Br"), Some(Element::Br));
        assert_eq!(Element::from_symbol("Xx"), None);
    }

    #[test]
    fn test_default_valences() {
        assert_eq!(Element::C.default_valences(), &[4]);
        assert_eq!(Element::N.default_valences(), &[3, 5]);
        assert_eq!(Element::S.default_valences(), &[2, 4, 6]);
    }

    #[test]
    fn test_is_organic_subset() {
        assert!(Element::C.is_organic_subset());
        assert!(!Element::Fe.is_organic_subset());
    }

    #[test]
    fn test_atomic_number() {
        assert_eq!(Element::C.atomic_number(), 6);
        assert_eq!(Element::N.atomic_number(), 7);
        assert_eq!(Element::Fe.atomic_number(), 26);
    }
}
