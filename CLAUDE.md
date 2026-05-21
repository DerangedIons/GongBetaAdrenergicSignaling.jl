# Julia Best Practices

## Julia Style & Naming
- `snake_case` for functions and variables, `CamelCase` for types and modules, `SCREAMING_SNAKE_CASE` for constants
- `!` suffix for mutating functions (e.g., `normalize!`)
- No type piracy — don't add methods to types you don't own for functions you don't own

## Performance
- Type stability: never change a variable's type mid-function; use concrete types in struct fields
- `const` for global variables, or pass them as function arguments
- Use `@views` for array slices to avoid allocations
- Pre-allocate outputs; prefer in-place operations (`mul!`, `ldiv!`, `@.`)
- Column-major iteration order — inner loop over first index
- No `Vector{Any}` or abstract element types in containers
- Performance-critical code must live inside functions, never at global scope

## Type System & Dispatch
- Prefer multiple dispatch over if/else type-checking
- Don't over-constrain argument types — use duck typing (e.g., `f(x)` not `f(x::Float64)`) unless dispatch requires it
- Use parametric types for generic, performant structs
- Abstract types for API/dispatch hierarchies; concrete types for storage

## Error Handling
- Use specific exception types (`ArgumentError`, `DimensionMismatch`, `BoundsError`, etc.)
- `throw` for actual errors, not control flow

## Package & Module Conventions
- Standard layout: `src/`, `test/`, `docs/`
- Explicit `export` — only export the public API
- Triple-quote docstrings (`"""..."""`) above functions

## Testing
- Investigate the project's test setup first (`Test`, `TestItemRunner`, etc.)
- Use `@test_throws` for error cases
- Each `@testitem` is self-contained with its own `using` statements (TestItems.jl requirement)
