use crate::mol::atom::Element;
use crate::mol::bond::BondType;
use thiserror::Error;

#[derive(Debug, Clone, PartialEq)]
pub enum Token {
    /// Organic subset atom (B, C, N, O, P, S, F, Cl, Br, I; lowercase = aromatic).
    Atom {
        element: Element,
        aromatic: bool,
    },
    /// Bracket atom: `[NH4+]`, `[13C@@H]`, `[Fe+2]`
    BracketAtom {
        isotope: Option<u16>,
        element: Element,
        aromatic: bool,
        h_count: Option<u8>,
        charge: i8,
        map_num: Option<u16>,
        chirality: Option<Chirality>,
    },
    Bond(BondType),
    OpenBranch,
    CloseBranch,
    RingClosure(u8),
    Dot,
    BondStereo(BondStereo),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Chirality {
    Anticlockwise, // @
    Clockwise,     // @@
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BondStereo {
    Up,   // /
    Down, // \
}

#[derive(Error, Debug)]
pub enum LexerError {
    #[error("unexpected character '{0}' at position {1}")]
    UnexpectedChar(char, usize),
    #[error("unclosed bracket atom starting at position {0}")]
    UnclosedBracket(usize),
    #[error("invalid element symbol in bracket at position {0}")]
    InvalidElement(usize),
    #[error("invalid charge specification at position {0}")]
    InvalidCharge(usize),
    #[error("invalid ring closure number at position {0}")]
    InvalidRingClosure(usize),
}

/// Tokenize a SMILES string into a vector of tokens.
pub fn tokenize(smiles: &str) -> Result<Vec<Token>, LexerError> {
    let bytes = smiles.as_bytes();
    let mut tokens = Vec::with_capacity(smiles.len());
    let mut i = 0;

    while i < bytes.len() {
        match bytes[i] {
            b'(' => {
                tokens.push(Token::OpenBranch);
                i += 1;
            }
            b')' => {
                tokens.push(Token::CloseBranch);
                i += 1;
            }
            b'.' => {
                tokens.push(Token::Dot);
                i += 1;
            }
            b'-' => {
                tokens.push(Token::Bond(BondType::Single));
                i += 1;
            }
            b'=' => {
                tokens.push(Token::Bond(BondType::Double));
                i += 1;
            }
            b'#' => {
                tokens.push(Token::Bond(BondType::Triple));
                i += 1;
            }
            b':' => {
                tokens.push(Token::Bond(BondType::Aromatic));
                i += 1;
            }
            b'/' => {
                tokens.push(Token::BondStereo(BondStereo::Up));
                i += 1;
            }
            b'\\' => {
                tokens.push(Token::BondStereo(BondStereo::Down));
                i += 1;
            }

            b'0'..=b'9' => {
                tokens.push(Token::RingClosure(bytes[i] - b'0'));
                i += 1;
            }
            b'%' => {
                if i + 2 < bytes.len()
                    && bytes[i + 1].is_ascii_digit()
                    && bytes[i + 2].is_ascii_digit()
                {
                    let num = (bytes[i + 1] - b'0') * 10 + (bytes[i + 2] - b'0');
                    tokens.push(Token::RingClosure(num));
                    i += 3;
                } else {
                    return Err(LexerError::InvalidRingClosure(i));
                }
            }

            b'[' => {
                let (token, consumed) = parse_bracket_atom(bytes, i)?;
                tokens.push(token);
                i += consumed;
            }

            b'B' => {
                if i + 1 < bytes.len() && bytes[i + 1] == b'r' {
                    tokens.push(Token::Atom {
                        element: Element::Br,
                        aromatic: false,
                    });
                    i += 2;
                } else {
                    tokens.push(Token::Atom {
                        element: Element::B,
                        aromatic: false,
                    });
                    i += 1;
                }
            }
            b'C' => {
                if i + 1 < bytes.len() && bytes[i + 1] == b'l' {
                    tokens.push(Token::Atom {
                        element: Element::Cl,
                        aromatic: false,
                    });
                    i += 2;
                } else {
                    tokens.push(Token::Atom {
                        element: Element::C,
                        aromatic: false,
                    });
                    i += 1;
                }
            }
            b'N' => {
                tokens.push(Token::Atom {
                    element: Element::N,
                    aromatic: false,
                });
                i += 1;
            }
            b'O' => {
                tokens.push(Token::Atom {
                    element: Element::O,
                    aromatic: false,
                });
                i += 1;
            }
            b'P' => {
                tokens.push(Token::Atom {
                    element: Element::P,
                    aromatic: false,
                });
                i += 1;
            }
            b'S' => {
                if i + 1 < bytes.len() && bytes[i + 1] == b'i' {
                    tokens.push(Token::Atom {
                        element: Element::Si,
                        aromatic: false,
                    });
                    i += 2;
                } else {
                    tokens.push(Token::Atom {
                        element: Element::S,
                        aromatic: false,
                    });
                    i += 1;
                }
            }
            b'F' => {
                tokens.push(Token::Atom {
                    element: Element::F,
                    aromatic: false,
                });
                i += 1;
            }
            b'I' => {
                tokens.push(Token::Atom {
                    element: Element::I,
                    aromatic: false,
                });
                i += 1;
            }

            // Aromatic organic subset atoms
            b'c' => {
                tokens.push(Token::Atom {
                    element: Element::C,
                    aromatic: true,
                });
                i += 1;
            }
            b'n' => {
                tokens.push(Token::Atom {
                    element: Element::N,
                    aromatic: true,
                });
                i += 1;
            }
            b'o' => {
                tokens.push(Token::Atom {
                    element: Element::O,
                    aromatic: true,
                });
                i += 1;
            }
            b's' => {
                if i + 1 < bytes.len() && bytes[i + 1] == b'e' {
                    tokens.push(Token::Atom {
                        element: Element::Se,
                        aromatic: true,
                    });
                    i += 2;
                } else {
                    tokens.push(Token::Atom {
                        element: Element::S,
                        aromatic: true,
                    });
                    i += 1;
                }
            }
            b'p' => {
                tokens.push(Token::Atom {
                    element: Element::P,
                    aromatic: true,
                });
                i += 1;
            }
            b'b' => {
                tokens.push(Token::Atom {
                    element: Element::B,
                    aromatic: true,
                });
                i += 1;
            }

            other => {
                return Err(LexerError::UnexpectedChar(other as char, i));
            }
        }
    }

    Ok(tokens)
}

/// Parse a bracket atom starting at `start` (points to '[').
/// Returns (Token, bytes consumed including '[' and ']').
fn parse_bracket_atom(bytes: &[u8], start: usize) -> Result<(Token, usize), LexerError> {
    let mut i = start + 1; // skip '['

    // Find closing ']' to bound parsing
    let end = bytes[i..]
        .iter()
        .position(|&b| b == b']')
        .ok_or(LexerError::UnclosedBracket(start))?
        + i;

    // Parse isotope (leading digits)
    let mut isotope: Option<u16> = None;
    if i < end && bytes[i].is_ascii_digit() {
        let mut num: u16 = 0;
        while i < end && bytes[i].is_ascii_digit() {
            num = num * 10 + (bytes[i] - b'0') as u16;
            i += 1;
        }
        isotope = Some(num);
    }

    // Parse element symbol
    if i >= end {
        return Err(LexerError::InvalidElement(i));
    }

    // Try two-char element first, then one-char, then aromatic lowercase
    let (element, aromatic, consumed) = if bytes[i].is_ascii_uppercase() {
        if i + 1 < end && bytes[i + 1].is_ascii_lowercase() && bytes[i + 1] != b'H' {
            // two-char element like Fe, Cl, Br, etc.
            let sym = std::str::from_utf8(&bytes[i..i + 2]).unwrap_or("");
            if let Some(el) = Element::from_symbol(sym) {
                (el, false, 2)
            } else {
                // Fall back to single char
                let sym1 = std::str::from_utf8(&bytes[i..i + 1]).unwrap_or("");
                let el = Element::from_symbol(sym1).ok_or(LexerError::InvalidElement(i))?;
                (el, false, 1)
            }
        } else {
            let sym = std::str::from_utf8(&bytes[i..i + 1]).unwrap_or("");
            let el = Element::from_symbol(sym).ok_or(LexerError::InvalidElement(i))?;
            (el, false, 1)
        }
    } else if bytes[i].is_ascii_lowercase() {
        // aromatic atom in bracket: c, n, o, s, p, se, te, etc.
        // Try two-char element first (e.g. [se] -> Se, [te] -> Te)
        if i + 1 < end && bytes[i + 1].is_ascii_lowercase() {
            let two = [bytes[i].to_ascii_uppercase(), bytes[i + 1]];
            let sym2 = std::str::from_utf8(&two).unwrap_or("");
            if let Some(el) = Element::from_symbol(sym2) {
                (el, true, 2)
            } else {
                let upper = bytes[i].to_ascii_uppercase();
                let sym = std::str::from_utf8(std::slice::from_ref(&upper)).unwrap_or("");
                let el = Element::from_symbol(sym).ok_or(LexerError::InvalidElement(i))?;
                (el, true, 1)
            }
        } else {
            let upper = bytes[i].to_ascii_uppercase();
            let sym = std::str::from_utf8(std::slice::from_ref(&upper)).unwrap_or("");
            let el = Element::from_symbol(sym).ok_or(LexerError::InvalidElement(i))?;
            (el, true, 1)
        }
    } else {
        return Err(LexerError::InvalidElement(i));
    };
    i += consumed;

    // Parse chirality (@ or @@)
    let chirality = if i < end && bytes[i] == b'@' {
        i += 1;
        if i < end && bytes[i] == b'@' {
            i += 1;
            Some(Chirality::Clockwise)
        } else {
            Some(Chirality::Anticlockwise)
        }
    } else {
        None
    };

    // Parse H count
    let h_count = if i < end && bytes[i] == b'H' {
        i += 1;
        if i < end && bytes[i].is_ascii_digit() {
            let h = bytes[i] - b'0';
            i += 1;
            Some(h)
        } else {
            Some(1u8)
        }
    } else {
        None
    };

    // Parse charge
    let charge: i8 = if i < end && (bytes[i] == b'+' || bytes[i] == b'-') {
        let sign: i8 = if bytes[i] == b'+' { 1 } else { -1 };
        i += 1;
        if i < end && bytes[i].is_ascii_digit() {
            // Numeric charge: +2, -3, etc.
            let mut mag: i8 = 0;
            while i < end && bytes[i].is_ascii_digit() {
                mag = mag * 10 + (bytes[i] - b'0') as i8;
                i += 1;
            }
            sign * mag
        } else if i < end && bytes[i] == bytes[i - 1] {
            // Repeated: ++ or --
            let mut mag: i8 = 1;
            while i < end && (bytes[i] == b'+' || bytes[i] == b'-') {
                mag += 1;
                i += 1;
            }
            sign * mag
        } else {
            sign
        }
    } else {
        0
    };

    // Parse atom map number (:N)
    let map_num = if i < end && bytes[i] == b':' {
        i += 1;
        let mut num: u16 = 0;
        while i < end && bytes[i].is_ascii_digit() {
            num = num * 10 + (bytes[i] - b'0') as u16;
            i += 1;
        }
        Some(num)
    } else {
        None
    };

    // i should now be at end (the ']')
    let consumed_total = end - start + 1; // includes '[' and ']'

    Ok((
        Token::BracketAtom {
            isotope,
            element,
            aromatic,
            h_count,
            charge,
            map_num,
            chirality,
        },
        consumed_total,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_organic() {
        let tokens = tokenize("CCO").unwrap();
        assert_eq!(tokens.len(), 3);
        assert_eq!(
            tokens[0],
            Token::Atom {
                element: Element::C,
                aromatic: false
            }
        );
        assert_eq!(
            tokens[1],
            Token::Atom {
                element: Element::C,
                aromatic: false
            }
        );
        assert_eq!(
            tokens[2],
            Token::Atom {
                element: Element::O,
                aromatic: false
            }
        );
    }

    #[test]
    fn test_two_char_elements() {
        let tokens = tokenize("ClBr").unwrap();
        assert_eq!(tokens.len(), 2);
        assert_eq!(
            tokens[0],
            Token::Atom {
                element: Element::Cl,
                aromatic: false
            }
        );
        assert_eq!(
            tokens[1],
            Token::Atom {
                element: Element::Br,
                aromatic: false
            }
        );
    }

    #[test]
    fn test_benzene_aromatic() {
        let tokens = tokenize("c1ccccc1").unwrap();
        assert_eq!(tokens.len(), 8);
        assert_eq!(
            tokens[0],
            Token::Atom {
                element: Element::C,
                aromatic: true
            }
        );
        assert_eq!(tokens[1], Token::RingClosure(1));
    }

    #[test]
    fn test_branch() {
        let tokens = tokenize("CC(=O)O").unwrap();
        assert_eq!(tokens.len(), 7);
        assert_eq!(tokens[2], Token::OpenBranch);
        assert_eq!(tokens[3], Token::Bond(BondType::Double));
        assert_eq!(tokens[5], Token::CloseBranch);
    }

    #[test]
    fn test_bracket_atom_simple() {
        let tokens = tokenize("[NH4+]").unwrap();
        assert_eq!(tokens.len(), 1);
        match &tokens[0] {
            Token::BracketAtom {
                element,
                h_count,
                charge,
                ..
            } => {
                assert_eq!(*element, Element::N);
                assert_eq!(*h_count, Some(4));
                assert_eq!(*charge, 1);
            }
            _ => panic!("expected BracketAtom"),
        }
    }

    #[test]
    fn test_bracket_atom_isotope() {
        let tokens = tokenize("[13C]").unwrap();
        match &tokens[0] {
            Token::BracketAtom {
                isotope, element, ..
            } => {
                assert_eq!(*isotope, Some(13));
                assert_eq!(*element, Element::C);
            }
            _ => panic!("expected BracketAtom"),
        }
    }

    #[test]
    fn test_bracket_atom_charge_negative() {
        let tokens = tokenize("[O-]").unwrap();
        match &tokens[0] {
            Token::BracketAtom {
                element, charge, ..
            } => {
                assert_eq!(*element, Element::O);
                assert_eq!(*charge, -1);
            }
            _ => panic!("expected BracketAtom"),
        }
    }

    #[test]
    fn test_bracket_atom_chirality() {
        let tokens = tokenize("[C@@H]").unwrap();
        match &tokens[0] {
            Token::BracketAtom {
                chirality, h_count, ..
            } => {
                assert_eq!(*chirality, Some(Chirality::Clockwise));
                assert_eq!(*h_count, Some(1));
            }
            _ => panic!("expected BracketAtom"),
        }
    }

    #[test]
    fn test_ring_closure_two_digit() {
        let tokens = tokenize("C%10CC%10").unwrap();
        assert_eq!(tokens[1], Token::RingClosure(10));
        assert_eq!(tokens[4], Token::RingClosure(10));
    }

    #[test]
    fn test_bond_stereo() {
        // "F/C=C/F" → [F, /(stereo), C, =(bond), C, /(stereo), F]
        let tokens = tokenize("F/C=C/F").unwrap();
        assert_eq!(tokens.len(), 7);
        assert_eq!(tokens[1], Token::BondStereo(BondStereo::Up));
        assert_eq!(tokens[5], Token::BondStereo(BondStereo::Up));
    }

    #[test]
    fn test_disconnected() {
        let tokens = tokenize("CC.OO").unwrap();
        assert_eq!(tokens[2], Token::Dot);
    }

    #[test]
    fn test_complex_drug() {
        let tokens = tokenize("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        assert!(tokens.len() > 15);
    }

    #[test]
    fn test_charged_iron() {
        let tokens = tokenize("[Fe+2]").unwrap();
        match &tokens[0] {
            Token::BracketAtom {
                element, charge, ..
            } => {
                assert_eq!(*element, Element::Fe);
                assert_eq!(*charge, 2);
            }
            _ => panic!("expected BracketAtom"),
        }
    }
}
