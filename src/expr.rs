use std::collections::HashMap;
use std::fmt::{self, Display, Formatter};
use std::str::FromStr;

#[derive(Clone, Debug, PartialEq)]
pub enum Expr {
    Constant(f64),
    Concentration(usize),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Div(Box<Expr>, Box<Expr>),
    Pow(Box<Expr>, Box<Expr>),
    Exp(Box<Expr>),
}

impl Expr {
    pub fn eval(&self, species: &[isize]) -> f64 {
        match self {
            Expr::Constant(c) => *c,
            Expr::Concentration(i) => *unsafe { species.get_unchecked(*i) } as f64,
            Expr::Add(a, b) => a.eval(species) + b.eval(species),
            Expr::Sub(a, b) => a.eval(species) - b.eval(species),
            Expr::Mul(a, b) => a.eval(species) * b.eval(species),
            Expr::Div(a, b) => a.eval(species) / b.eval(species),
            Expr::Pow(a, b) => a.eval(species).powf(b.eval(species)),
            Expr::Exp(a) => a.eval(species).exp(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum PExpr {
    Constant(f64),
    Variable(String),
    Add(Box<PExpr>, Box<PExpr>),
    Sub(Box<PExpr>, Box<PExpr>),
    Mul(Box<PExpr>, Box<PExpr>),
    Div(Box<PExpr>, Box<PExpr>),
    Pow(Box<PExpr>, Box<PExpr>),
    Exp(Box<PExpr>),
}

impl PExpr {
    pub fn to_expr(
        &self,
        species: &HashMap<String, usize>,
        params: &HashMap<String, f64>,
    ) -> Result<Expr, String> {
        let expr = match self {
            PExpr::Constant(c) => Expr::Constant(*c),
            PExpr::Variable(s) => {
                if let Some(i) = species.get(s) {
                    Expr::Concentration(*i)
                } else {
                    let v = params
                        .get(s)
                        .ok_or(format!("Parameter {s} should have a value").to_string())?;
                    Expr::Constant(*v)
                }
            }
            PExpr::Add(a, b) => Expr::Add(
                Box::new(a.to_expr(species, params)?),
                Box::new(b.to_expr(species, params)?),
            ),
            PExpr::Sub(a, b) => Expr::Sub(
                Box::new(a.to_expr(species, params)?),
                Box::new(b.to_expr(species, params)?),
            ),
            PExpr::Mul(a, b) => Expr::Mul(
                Box::new(a.to_expr(species, params)?),
                Box::new(b.to_expr(species, params)?),
            ),
            PExpr::Div(a, b) => Expr::Div(
                Box::new(a.to_expr(species, params)?),
                Box::new(b.to_expr(species, params)?),
            ),
            PExpr::Pow(a, b) => Expr::Pow(
                Box::new(a.to_expr(species, params)?),
                Box::new(b.to_expr(species, params)?),
            ),
            PExpr::Exp(a) => Expr::Exp(Box::new(a.to_expr(species, params)?)),
        };
        Ok(expr)
    }
}

impl Display for PExpr {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        match self {
            PExpr::Constant(c) => write!(f, "{c}"),
            PExpr::Variable(c) => write!(f, "{c}"),
            PExpr::Add(a, b) => write!(f, "({a} + {b})"),
            PExpr::Sub(a, b) => write!(f, "({a} - {b})"),
            PExpr::Mul(a, b) => write!(f, "({a} * {b})"),
            PExpr::Div(a, b) => write!(f, "({a} / {b})"),
            PExpr::Pow(a, b) => write!(f, "({a} ^ {b})"),
            PExpr::Exp(a) => write!(f, "exp({a})"),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct RateParseError;

impl FromStr for PExpr {
    type Err = RateParseError;

    fn from_str(mut s: &str) -> Result<Self, Self::Err> {
        parsing::nomexpr(&mut s).map_err(|_| RateParseError {})
    }
}

mod parsing {
    use crate::expr::PExpr;
    use winnow::ascii::{float, space0};
    use winnow::combinator::{alt, delimited, preceded, separated_foldl1, separated_pair};
    use winnow::error::{ContextError, ParseError};
    use winnow::prelude::*;
    use winnow::token::{one_of, take_while};
    use winnow::Result;

    fn constant(s: &mut &str) -> Result<PExpr> {
        float.map(PExpr::Constant).parse_next(s)
    }

    fn variable(s: &mut &str) -> Result<PExpr> {
        // from winnow recipes
        // https://docs.rs/winnow/latest/winnow/_topic/language/index.html
        (
            one_of(|c: char| c.is_alpha() || c == '_'),
            take_while(0.., |c: char| c.is_alphanum() || c == '_'),
        )
            .take()
            .map(str::to_string)
            .map(PExpr::Variable)
            .parse_next(s)
    }

    fn parentheses(s: &mut &str) -> Result<PExpr> {
        delimited("(", delimited(space0, expr, space0), ")").parse_next(s)
    }

    fn expr(s: &mut &str) -> Result<PExpr> {
        alt((add_sub, term)).parse_next(s)
    }

    fn term(s: &mut &str) -> Result<PExpr> {
        alt((mul_div, factor)).parse_next(s)
    }

    fn factor(s: &mut &str) -> Result<PExpr> {
        alt((pow, atom)).parse_next(s)
    }

    fn atom(s: &mut &str) -> Result<PExpr> {
        // variable must be before constant
        // otherwise it matches variables starting with inf or nan as a float
        alt((exp, variable, constant, parentheses)).parse_next(s)
    }

    fn add_sub(s: &mut &str) -> Result<PExpr> {
        separated_foldl1(
            term,
            delimited(space0, one_of(['+', '-']), space0),
            |a, op, b| match op {
                '+' => PExpr::Add(Box::new(a), Box::new(b)),
                '-' => PExpr::Sub(Box::new(a), Box::new(b)),
                _ => unreachable!(),
            },
        )
        .parse_next(s)
    }

    fn mul_div(s: &mut &str) -> Result<PExpr> {
        separated_foldl1(
            factor,
            delimited(space0, one_of(['*', '/']), space0),
            |a, op, b| match op {
                '*' => PExpr::Mul(Box::new(a), Box::new(b)),
                '/' => PExpr::Div(Box::new(a), Box::new(b)),
                _ => unreachable!(),
            },
        )
        .parse_next(s)
    }

    fn pow(s: &mut &str) -> Result<PExpr> {
        separated_pair(atom, (space0, '^', space0), atom)
            .map(|(a, b)| PExpr::Pow(Box::new(a), Box::new(b)))
            .parse_next(s)
    }

    fn exp(s: &mut &str) -> Result<PExpr> {
        preceded("exp", parentheses)
            .map(Box::new)
            .map(PExpr::Exp)
            .parse_next(s)
    }

    pub(crate) fn nomexpr<'a>(s: &mut &'a str) -> Result<PExpr, ParseError<&'a str, ContextError>> {
        expr.parse(s)
    }

    #[cfg(test)]
    mod tests {
        use crate::expr::{PExpr, RateParseError};

        #[test]
        fn test_constant() {
            let e: PExpr = "3".parse().unwrap();
            assert_eq!(format!("{e}"), "3");
            let e: PExpr = ".23".parse().unwrap();
            assert_eq!(format!("{e}"), "0.23");
            let e: PExpr = "1.23".parse().unwrap();
            assert_eq!(format!("{e}"), "1.23");
            let e: PExpr = "2e-3".parse().unwrap();
            assert_eq!(format!("{e}"), "0.002");
            let e: PExpr = "2.04e3".parse().unwrap();
            assert_eq!(format!("{e}"), "2040");
            let e: PExpr = "-2.04e3".parse().unwrap();
            assert_eq!(format!("{e}"), "-2040");
            let e: PExpr = "+2.04e3".parse().unwrap();
            assert_eq!(format!("{e}"), "2040");
        }

        #[test]
        fn test_concentration() {
            let e: PExpr = "a".parse().unwrap();
            assert_eq!(format!("{e}"), "a");
            let e: PExpr = "_".parse().unwrap();
            assert_eq!(format!("{e}"), "_");
            let e: PExpr = "exp".parse().unwrap();
            assert_eq!(format!("{e}"), "exp");
            let e: PExpr = "PhoP3".parse().unwrap();
            assert_eq!(format!("{e}"), "PhoP3");
            let e: PExpr = "Mono1_Mono2_".parse().unwrap();
            assert_eq!(format!("{e}"), "Mono1_Mono2_");
        }

        #[test]
        fn test_exp() {
            let e: PExpr = "exp(3)".parse().unwrap();
            assert_eq!(format!("{e}"), "exp(3)");
        }

        #[test]
        fn test_pow() {
            let e: PExpr = "3^4".parse().unwrap();
            assert_eq!(format!("{e}"), "(3 ^ 4)");
            let e: PExpr = "3 ^4".parse().unwrap();
            assert_eq!(format!("{e}"), "(3 ^ 4)");
            let e: PExpr = "3^ 4".parse().unwrap();
            assert_eq!(format!("{e}"), "(3 ^ 4)");
            let e: PExpr = "3 ^ 4".parse().unwrap();
            assert_eq!(format!("{e}"), "(3 ^ 4)");
        }

        #[test]
        fn test_mul_div() {
            let e: PExpr = "3*4".parse().unwrap();
            assert_eq!(format!("{e}"), "(3 * 4)");
            let e: PExpr = "3 *4".parse().unwrap();
            assert_eq!(format!("{e}"), "(3 * 4)");
            let e: PExpr = "3* 4".parse().unwrap();
            assert_eq!(format!("{e}"), "(3 * 4)");
            let e: PExpr = "3 * 4".parse().unwrap();
            assert_eq!(format!("{e}"), "(3 * 4)");
            let e: PExpr = "3*4/2".parse().unwrap();
            assert_eq!(format!("{e}"), "((3 * 4) / 2)");
            let e: PExpr = "3*4 /2".parse().unwrap();
            assert_eq!(format!("{e}"), "((3 * 4) / 2)");
            let e: PExpr = "3*4/ 2".parse().unwrap();
            assert_eq!(format!("{e}"), "((3 * 4) / 2)");
            let e: PExpr = "3*4 / 2".parse().unwrap();
            assert_eq!(format!("{e}"), "((3 * 4) / 2)");
        }

        #[test]
        fn test_add_sub() {
            let e: PExpr = "3+4".parse().unwrap();
            assert_eq!(format!("{e}"), "(3 + 4)");
            let e: PExpr = "3-4+1".parse().unwrap();
            assert_eq!(format!("{e}"), "((3 - 4) + 1)");
            let e: PExpr = "3-4-1".parse().unwrap();
            assert_eq!(format!("{e}"), "((3 - 4) - 1)");
        }

        #[test]
        fn test_expr() {
            let e: PExpr = "1.20 * A*B".parse().unwrap();
            assert_eq!(format!("{e}"), "((1.2 * A) * B)");
            let e: PExpr = "1.20*Sugar / (3.5+Sugar)".parse().unwrap();
            assert_eq!(format!("{e}"), "((1.2 * Sugar) / (3.5 + Sugar))");
        }

        #[test]
        fn test_fail() {
            let e: Result<PExpr, _> = "+".parse();
            assert_eq!(e, Result::Err(RateParseError));
            let e: Result<PExpr, _> = "1+".parse();
            assert_eq!(e, Result::Err(RateParseError));
        }

        #[test]
        fn test_float_names() {
            // inf
            let e: PExpr = "inf".parse().unwrap();
            assert_eq!(PExpr::Variable("inf".to_string()), e);
            // starts with inf
            let e: PExpr = "infect".parse().unwrap();
            assert_eq!(PExpr::Variable("infect".to_string()), e);
            // nan
            let e: PExpr = "nan".parse().unwrap();
            assert_eq!(PExpr::Variable("nan".to_string()), e);
            // starts with nan
            let e: PExpr = "nanny".parse().unwrap();
            assert_eq!(PExpr::Variable("nanny".to_string()), e);
            // E
            let e: PExpr = "E".parse().unwrap();
            assert_eq!(PExpr::Variable("E".to_string()), e);
            // e
            let e: PExpr = "e".parse().unwrap();
            assert_eq!(PExpr::Variable("e".to_string()), e);
            // exp
            let e: PExpr = "exp".parse().unwrap();
            assert_eq!(PExpr::Variable("exp".to_string()), e);
            // starts with exp
            let e: PExpr = "explicit".parse().unwrap();
            assert_eq!(PExpr::Variable("explicit".to_string()), e);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::expr::{Expr, PExpr};
    use std::collections::HashMap;

    #[test]
    fn test_conversion() {
        let pe: PExpr = "1.21 * C + B - A / D ^ E * (F + exp(D))".parse().unwrap();
        // fmt: off
        let e = Expr::Sub(
            Box::new(Expr::Add(
                Box::new(Expr::Mul(
                    Box::new(Expr::Constant(1.21)),
                    Box::new(Expr::Concentration(2)),
                )),
                Box::new(Expr::Concentration(1)),
            )),
            Box::new(Expr::Mul(
                Box::new(Expr::Div(
                    Box::new(Expr::Concentration(0)),
                    Box::new(Expr::Pow(
                        Box::new(Expr::Concentration(3)),
                        Box::new(Expr::Concentration(4)),
                    )),
                )),
                Box::new(Expr::Add(
                    Box::new(Expr::Concentration(5)),
                    Box::new(Expr::Exp(Box::new(Expr::Concentration(3)))),
                )),
            )),
        );
        // fmt: on
        let mut init = HashMap::new();
        let mut params = HashMap::new();
        init.insert("A".to_string(), 0);
        init.insert("B".to_string(), 1);
        init.insert("C".to_string(), 2);
        init.insert("D".to_string(), 3);
        init.insert("E".to_string(), 4);
        init.insert("F".to_string(), 5);
        assert_eq!(pe.to_expr(&init, &params), Ok(e));
    }

    #[test]
    fn test_eval() {
        let pe: PExpr = "1.21 * C + B - A / D ^ E * (F + exp(D))".parse().unwrap();
        let mut init = HashMap::new();
        let mut params = HashMap::new();
        init.insert("A".to_string(), 0);
        init.insert("B".to_string(), 1);
        init.insert("C".to_string(), 2);
        init.insert("D".to_string(), 3);
        init.insert("E".to_string(), 4);
        init.insert("F".to_string(), 5);
        let mut species = vec![2, 3, 5, 7, 11, 13];
        assert_eq!(
            pe.to_expr(&init, &params).unwrap().eval(&species),
            9.049998877643098
        );
    }
}
