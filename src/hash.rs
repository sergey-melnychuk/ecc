use crate::{hex, Int};

pub fn hash(input: &[u8]) -> Int {
    hex(&kangarootwelve_xkcp::hash(input)
        .to_hex()
        .chars()
        .collect::<String>())
}

pub fn hash_n(inputs: &[&[u8]]) -> Int {
    let mut hasher = kangarootwelve_xkcp::Hasher::new();
    for input in inputs {
        hasher.update(input);
    }
    hex(&hasher.finalize().to_hex().chars().collect::<String>())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_hash() {
        let int = crate::hash::hash(b"it works");
        let exp = crate::dec(
            "114718067966976356955517098759415765868171761422675185365940576053772770837809",
        );
        assert_eq!(int, exp);
    }
}
