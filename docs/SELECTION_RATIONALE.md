# Selection rationale

The package is deliberately limited to functions needed to follow a reported
result or a sensitivity analysis that materially constrains interpretation.

Excluded material comprises:

- figure-generation and presentation code;
- server-specific runners and path-discovery scripts;
- generated outputs, logs, and intermediate workspaces;
- participant identifiers and identifier mappings;
- generic inventory-only audits whose numerical outputs are already supplied
  in the Source Data;
- separate early beta and theta-alpha edge functions superseded by the unified
  630-endpoint implementations;
- beta-subband orchestration functions that duplicate the public source,
  event-detection, and inference functions;
- internal release checklists and development notes.

The cropped-envelope edge function is retained because the manuscript reports
its participant-level agreement with the mirror-padded implementation. The
mirror-padded edge, spatial-leakage, and orthogonalization functions are all
retained because together they define the manuscript's limits on anatomical
interpretation.

