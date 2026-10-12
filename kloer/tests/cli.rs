use std::process::Command;

#[test]
fn exploratory_warning_is_shown_with_help_and_usage_errors() {
    for args in [
        vec!["--help"],
        vec!["compress", "--help"],
        vec!["compress-list", "--help"],
        vec!["decompress", "--help"],
        vec!["decompress-list", "--help"],
        vec![],
        vec!["--unknown-option"],
        vec!["compress"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_kloer"))
            .args(&args)
            .output()
            .unwrap();
        let text = format!(
            "{}{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            text.contains("WARNING: KLOER is an exploratory project."),
            "{args:?}: {text}"
        );
        assert!(text.contains("Usage:"), "{args:?}: {text}");
        assert_eq!(output.status.success(), args.contains(&"--help"));
    }
}
