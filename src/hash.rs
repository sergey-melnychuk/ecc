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
    use super::*;
    #[test]
    fn test_hash() {
        let expected = "fda02020fb19df5f442aa940a4c4c4e89bcc7f2bf1ea586f408973844fd78531";
        assert_eq!(hash(b"it works"), hex(expected));
    }
}
