# Coding guidelines

## General
- Never use MATLAB classes.
- Use functions and structs only.
- Keep functions short.
- Do not change public interfaces unless necessary.
- Preserve backwards compatibility.

## Style
- Use lower_case function names.
- Comment all public functions.
- Keep line length below 100 characters.

## Philosophy
- Prioritize readability over clever code.
- Avoid unnecessary abstractions.
- Make minimal changes to existing code.

## Bayesian modelling
- Never modify the inference algorithm without discussion.
- Hyperparameters should be estimated by evidence optimization.
