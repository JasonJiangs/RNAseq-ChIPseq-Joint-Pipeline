# Framework Refactor

To make the framework more extensible, we need to refactor the framework with design patterns.

## Refactor Steps
### RNA-seq and ChIP-seq Tools
1. For each tool, we need to create a class that implements the `ITool` interface.
2. For each tool, we need to create a class that implements the `IToolFactory` interface.
3. For each parameter of a tool, we need to create a class that implements the `IParameter` interface.
### Machine Learning Tools
1. For each tool, we need to create a class that implements the `IMLTool` interface.
2. For each tool, we need to create a class that implements the `IMLToolFactory` interface.
3. For each parameter of a tool, we need to create a class that implements the `IMLParameter` interface.
### Tool Manager
1. Create a class that implements the `IToolManager` interface.
2. Create a class that implements the `IMLToolManager` interface.


## Design Patterns
https://sourcemaking.com/design_patterns
