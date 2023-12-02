use std::collections::HashMap;
use std::fmt::{self, Display, Formatter};
use std::str::FromStr;

#[derive(Clone, Debug)]
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

#[derive(Clone, Debug)]
pub(crate) enum NomExpr {
    Constant(f64),
    Concentration(String),
    Add(Box<NomExpr>, Box<NomExpr>),
    Sub(Box<NomExpr>, Box<NomExpr>),
    Mul(Box<NomExpr>, Box<NomExpr>),
    Div(Box<NomExpr>, Box<NomExpr>),
    Pow(Box<NomExpr>, Box<NomExpr>),
    Exp(Box<NomExpr>),
}

impl NomExpr {
    pub fn to_expr(&self, species: &HashMap<String, usize>) -> Expr {
        match self {
            NomExpr::Constant(c) => Expr::Constant(*c),
            NomExpr::Concentration(s) => Expr::Concentration(species[s]),
            NomExpr::Add(a, b) => {
                Expr::Add(Box::new(a.to_expr(species)), Box::new(b.to_expr(species)))
            }
            NomExpr::Sub(a, b) => {
                Expr::Sub(Box::new(a.to_expr(species)), Box::new(b.to_expr(species)))
            }
            NomExpr::Mul(a, b) => {
                Expr::Mul(Box::new(a.to_expr(species)), Box::new(b.to_expr(species)))
            }
            NomExpr::Div(a, b) => {
                Expr::Div(Box::new(a.to_expr(species)), Box::new(b.to_expr(species)))
            }
            NomExpr::Pow(a, b) => {
                Expr::Pow(Box::new(a.to_expr(species)), Box::new(b.to_expr(species)))
            }
            NomExpr::Exp(a) => Expr::Exp(Box::new(a.to_expr(species))),
        }
    }
}

impl Display for NomExpr {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        match self {
            NomExpr::Constant(c) => write!(f, "{c}"),
            NomExpr::Concentration(c) => write!(f, "{c}"),
            NomExpr::Add(a, b) => write!(f, "({a} + {b})"),
            NomExpr::Sub(a, b) => write!(f, "({a} - {b})"),
            NomExpr::Mul(a, b) => write!(f, "({a} * {b})"),
            NomExpr::Div(a, b) => write!(f, "({a} / {b})"),
            NomExpr::Pow(a, b) => write!(f, "({a} ^ {b})"),
            NomExpr::Exp(a) => write!(f, "exp({a})"),
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct RateParseError;

impl FromStr for NomExpr {
    type Err = RateParseError;

    fn from_str(mut s: &str) -> Result<Self, Self::Err> {
        parsing::nomexpr(&mut s).map_err(|_| RateParseError {})
    }
}

mod parsing {
    use super::NomExpr;
    use winnow::ascii::{alpha1, alphanumeric1, float, space0};
    use winnow::combinator::{alt, delimited, fold_repeat, preceded, repeat};
    use winnow::prelude::*;
    use winnow::token::one_of;
    use winnow::Parser;

    fn constant(s: &mut &str) -> PResult<NomExpr> {
        float.map(|c| NomExpr::Constant(c)).parse_next(s)
    }

    fn concentration(s: &mut &str) -> PResult<NomExpr> {
        // from nom recipes
        (
            alt((alpha1, "_")),
            repeat::<_, _, Vec<_>, _, _>(0.., alt((alphanumeric1, "_"))),
        )
            .recognize()
            .map(|c: &str| NomExpr::Concentration(c.to_string()))
            .parse_next(s)
    }

    fn parentheses(s: &mut &str) -> PResult<NomExpr> {
        delimited("(", delimited(space0, expr, space0), ")").parse_next(s)
    }

    fn expr(s: &mut &str) -> PResult<NomExpr> {
        alt((add_sub, term)).parse_next(s)
    }

    fn term(s: &mut &str) -> PResult<NomExpr> {
        alt((mul_div, factor)).parse_next(s)
    }

    fn factor(s: &mut &str) -> PResult<NomExpr> {
        alt((pow, atom)).parse_next(s)
    }

    fn atom(s: &mut &str) -> PResult<NomExpr> {
        alt((constant, exp, concentration, parentheses)).parse_next(s)
    }

    fn add_sub(s: &mut &str) -> PResult<NomExpr> {
        let first = term(s)?;
        let rest = fold_repeat(
            0..,
            (delimited(space0, one_of(['+', '-']), space0), term),
            || first.clone(),
            |acc: NomExpr, (op, f)| match op {
                '+' => NomExpr::Add(Box::new(acc), Box::new(f)),
                '-' => NomExpr::Sub(Box::new(acc), Box::new(f)),
                _ => unreachable!(),
            },
        )
        .parse_next(s)?;
        Ok(rest)
    }

    fn mul_div(s: &mut &str) -> PResult<NomExpr> {
        let first = factor(s)?;
        let rest = fold_repeat(
            0..,
            (delimited(space0, one_of(['*', '/']), space0), factor),
            || first.clone(),
            |acc: NomExpr, (op, f)| match op {
                '*' => NomExpr::Mul(Box::new(acc), Box::new(f)),
                '/' => NomExpr::Div(Box::new(acc), Box::new(f)),
                _ => unreachable!(),
            },
        )
        .parse_next(s)?;
        Ok(rest)
    }

    fn pow(s: &mut &str) -> PResult<NomExpr> {
        let a = atom(s)?;
        let _ = delimited(space0, "^", space0).parse_next(s)?;
        let b = atom(s)?;
        Ok(NomExpr::Pow(Box::new(a), Box::new(b)))
    }

    fn exp(s: &mut &str) -> PResult<NomExpr> {
        preceded("exp", parentheses)
            .map(|e| NomExpr::Exp(Box::new(e)))
            .parse_next(s)
    }

    pub(crate) fn nomexpr(s: &mut &str) -> PResult<NomExpr> {
        expr.parse_next(s)
    }

    #[test]
    fn test_constant() {
        let e: NomExpr = "3".parse().unwrap();
        assert_eq!(format!("{e}"), "3");
        let e: NomExpr = ".23".parse().unwrap();
        assert_eq!(format!("{e}"), "0.23");
        let e: NomExpr = "1.23".parse().unwrap();
        assert_eq!(format!("{e}"), "1.23");
        let e: NomExpr = "2e-3".parse().unwrap();
        assert_eq!(format!("{e}"), "0.002");
        let e: NomExpr = "2.04e3".parse().unwrap();
        assert_eq!(format!("{e}"), "2040");
        let e: NomExpr = "-2.04e3".parse().unwrap();
        assert_eq!(format!("{e}"), "-2040");
    }

    #[test]
    fn test_concentration() {
        let e: NomExpr = "a".parse().unwrap();
        assert_eq!(format!("{e}"), "a");
        let e: NomExpr = "_".parse().unwrap();
        assert_eq!(format!("{e}"), "_");
        let e: NomExpr = "exp".parse().unwrap();
        assert_eq!(format!("{e}"), "exp");
        let e: NomExpr = "PhoP3".parse().unwrap();
        assert_eq!(format!("{e}"), "PhoP3");
        let e: NomExpr = "Mono1_Mono2_".parse().unwrap();
        assert_eq!(format!("{e}"), "Mono1_Mono2_");
    }

    #[test]
    fn test_exp() {
        let e: NomExpr = "exp(3)".parse().unwrap();
        assert_eq!(format!("{e}"), "exp(3)");
    }

    #[test]
    fn test_pow() {
        let e: NomExpr = "3^4".parse().unwrap();
        assert_eq!(format!("{e}"), "(3 ^ 4)");
        let e: NomExpr = "3 ^4".parse().unwrap();
        assert_eq!(format!("{e}"), "(3 ^ 4)");
        let e: NomExpr = "3^ 4".parse().unwrap();
        assert_eq!(format!("{e}"), "(3 ^ 4)");
        let e: NomExpr = "3 ^ 4".parse().unwrap();
        assert_eq!(format!("{e}"), "(3 ^ 4)");
    }

    #[test]
    fn test_mul_div() {
        let e: NomExpr = "3*4".parse().unwrap();
        assert_eq!(format!("{e}"), "(3 * 4)");
        let e: NomExpr = "3 *4".parse().unwrap();
        assert_eq!(format!("{e}"), "(3 * 4)");
        let e: NomExpr = "3* 4".parse().unwrap();
        assert_eq!(format!("{e}"), "(3 * 4)");
        let e: NomExpr = "3 * 4".parse().unwrap();
        assert_eq!(format!("{e}"), "(3 * 4)");
        let e: NomExpr = "3*4/2".parse().unwrap();
        assert_eq!(format!("{e}"), "((3 * 4) / 2)");
        let e: NomExpr = "3*4 /2".parse().unwrap();
        assert_eq!(format!("{e}"), "((3 * 4) / 2)");
        let e: NomExpr = "3*4/ 2".parse().unwrap();
        assert_eq!(format!("{e}"), "((3 * 4) / 2)");
        let e: NomExpr = "3*4 / 2".parse().unwrap();
        assert_eq!(format!("{e}"), "((3 * 4) / 2)");
    }

    #[test]
    fn test_add_sub() {
        let e: NomExpr = "3+4".parse().unwrap();
        assert_eq!(format!("{e}"), "(3 + 4)");
        let e: NomExpr = "3-4+1".parse().unwrap();
        assert_eq!(format!("{e}"), "((3 - 4) + 1)");
        let e: NomExpr = "3-4-1".parse().unwrap();
        assert_eq!(format!("{e}"), "((3 - 4) - 1)");
    }

    #[test]
    fn test_expr() {
        let e: NomExpr = "1.20 * A*B".parse().unwrap();
        assert_eq!(format!("{e}"), "((1.2 * A) * B)");
        let e: NomExpr = "1.20*Sugar / (3.5+Sugar)".parse().unwrap();
        assert_eq!(format!("{e}"), "((1.2 * Sugar) / (3.5 + Sugar))");
    }
}
