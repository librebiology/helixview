use std::collections::BTreeSet;
use std::sync::Arc;

use iced::{Element, Task, Theme};
use iced::keyboard::{self, Key, key::Named, Modifiers};
use iced::widget::{text_editor, text_input};
use helixview_core::Alignment as SeqAlignment;
use helixview_core::Sequence as CoreSequence;
use helixview_core::{
    History as EditHistory,
    InsertGapColumn, DeleteColumn as DelColumn, MoveSequence, SetResidues,
    InsertGapInSeq, DeleteGapInSeq,
    EditFeature, AddFeature, DeleteFeature as DeleteFeatureCmd,
};

use helixview_analysis::PairwiseResult;
use helixview_external::ExternalAligner;
use crate::blast::{self, BlastHit};
use crate::entrez;
use crate::views::{alignment_view, welcome_view, identity_matrix_view, col_summary_view, dot_plot_view, seq_editor_view, color_editor_view, conservation_view};
use crate::prefs::Preferences;
use crate::color_table::ColorTable;

// ── Shaded export options ─────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct ShadedExportOptions {
    pub scheme:         ColorScheme,
    pub threshold:      f32,
    pub show_ruler:     bool,
    pub show_consensus: bool,
}

impl Default for ShadedExportOptions {
    fn default() -> Self {
        Self {
            scheme:         ColorScheme::Identity,
            threshold:      0.50,
            show_ruler:     true,
            show_consensus: true,
        }
    }
}

// ── Accessory app ─────────────────────────────────────────────────────────────

/// A configured external tool that can be run from within HelixView.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct AccessoryApp {
    /// Display name shown in the menu.
    pub name: String,
    /// Full shell command. Use `{input}` for the temp FASTA input file path,
    /// `{output}` for a temp output file path.
    pub command: String,
}

impl Default for AccessoryApp {
    fn default() -> Self {
        Self {
            name:    "My Tool".to_string(),
            command: "my_tool {input}".to_string(),
        }
    }
}

// ── BLAST state ───────────────────────────────────────────────────────────────

#[derive(Debug, Default, Clone)]
pub enum BlastState {
    #[default]
    Idle,
    /// Waiting for NCBI: RID is known, results not ready yet.
    Submitted(String),
    /// Results received and parsed.
    Complete(Vec<BlastHit>),
    /// Submission or polling failed.
    Failed(String),
}

// ── Tab record (background-tab snapshot) ─────────────────────────────────────

/// Saved state for a non-active tab.
pub struct TabRecord {
    pub document:          Option<Arc<SeqAlignment>>,
    pub view:              ViewState,
    pub history:           EditHistory,
    pub show_restr_map:    bool,
    pub protein_view_doc:  Option<Arc<SeqAlignment>>,
    pub protein_view_frame: usize,
}

impl TabRecord {
    fn empty() -> Self {
        Self {
            document:           None,
            view:               ViewState::default(),
            history:            EditHistory::new(100),
            show_restr_map:     false,
            protein_view_doc:   None,
            protein_view_frame: 0,
        }
    }
}

// ── Application state ─────────────────────────────────────────────────────────

pub struct HelixViewApp {
    /// Currently open alignment document.
    pub document: Option<Arc<SeqAlignment>>,
    /// Status message shown in the bottom bar.
    pub status: String,
    /// View state: scroll offsets, selection, etc.
    pub view: ViewState,
    /// Edit history for undo/redo.
    pub history: EditHistory,
    /// Most recent pairwise alignment result (shown in pairwise panel).
    pub pairwise_result: Option<PairwiseResult>,
    /// ORF finder results: (sequence name, list of ORFs).
    pub orf_results: Option<(String, Vec<helixview_analysis::Orf>)>,
    /// Whether the identity matrix view is open.
    pub show_identity_matrix: bool,
    /// Whether the positional column summary view is open.
    pub show_col_summary: bool,
    /// Sequence title right-click context menu state: Some(idx) when open.
    pub title_menu: Option<usize>,
    /// Text being typed into the rename field.
    pub rename_input: String,
    /// Whether the restriction map strip is visible.
    pub show_restr_map: bool,
    /// Whether the dot-plot view is open.
    pub show_dot_plot: bool,
    /// Dot-plot: index of the X-axis sequence.
    pub dot_plot_seq_a: usize,
    /// Dot-plot: index of the Y-axis sequence.
    pub dot_plot_seq_b: usize,
    /// Dot-plot: sliding window size (default 7).
    pub dot_plot_window: usize,
    /// Dot-plot: minimum matches in window to draw a dot (default 5).
    pub dot_plot_min_match: usize,
    /// Sequence raw-text editor: index of the sequence being edited (None = closed).
    pub seq_editor_idx: Option<usize>,
    /// Sequence raw-text editor: the mutable editor content (None when closed).
    pub seq_editor_content: Option<text_editor::Content>,
    /// Active residue color table (shared with the grid renderer).
    pub color_table: Arc<ColorTable>,
    /// Whether the color table editor view is open.
    pub show_color_editor: bool,
    /// Which residue is currently selected for editing: (byte, is_nuc).
    pub color_editor_selected: Option<(u8, bool)>,
    /// Whether the conservation search view is open.
    pub show_conservation: bool,
    /// Conservation search: identity threshold (0.0–1.0).
    pub conservation_threshold: f64,
    /// Conservation search: minimum consecutive conserved columns.
    pub conservation_min_width: usize,
    /// Six-frame translation results: (seq_name, nt_len, [frame_label, protein; 6]).
    pub six_frame_result: Option<(String, usize, [(i8, Vec<u8>); 6])>,
    /// Protein view: derived alignment shown instead of the DNA document (None = DNA mode).
    pub protein_view_doc: Option<Arc<SeqAlignment>>,
    /// Which reading frame is active for the protein view (0-based, 0–5).
    pub protein_view_frame: usize,
    /// BLAST state machine.
    pub blast: BlastState,
    /// All open tabs (background snapshots; active tab lives in the fields above).
    pub tabs: Vec<TabRecord>,
    /// Index of the currently active tab.
    pub active_tab: usize,
    /// Whether the composition analysis view is open.
    pub show_composition: bool,
    /// Whether the Tm calculator view is open.
    pub show_oligo_tm: bool,
    /// Phylogenetic tree viewer: currently loaded tree (None = closed).
    pub tree_data: Option<Arc<helixview_core::PhyloTree>>,
    /// Leaf names currently selected in the tree viewer.
    pub tree_selected: BTreeSet<String>,
    /// ABI trace viewer: currently loaded trace (None = trace view closed).
    pub trace_data: Option<Arc<helixview_formats::abi::AbiTrace>>,
    /// Horizontal scroll fraction for the trace viewer (0.0–1.0).
    pub trace_scroll: f32,
    /// Zoom level for the trace viewer (1.0 = all samples visible, up to 16.0).
    pub trace_zoom: f32,
    /// Quantification channel A (base letter, e.g. b'A').
    pub trace_quant_ch_a: u8,
    /// Quantification channel B (base letter, e.g. b'C').
    pub trace_quant_ch_b: u8,
    /// Loaded user preferences.
    pub prefs: Preferences,
    /// Whether the preferences dialog is open.
    pub show_prefs: bool,
    /// Whether the keyboard shortcut help view is open.
    pub show_help: bool,
    /// Whether the shaded graphic export dialog is open.
    pub show_shaded_export: bool,
    /// Options for the shaded graphic export.
    pub shaded_export_opts: ShadedExportOptions,
    /// Whether the taxonomy table view is open.
    pub show_taxonomy: bool,
    /// Whether the circular plasmid map view is open.
    pub show_plasmid: bool,
    /// Which sequence is displayed in the plasmid view.
    pub plasmid_seq_idx: usize,
    /// Whether feature arcs are drawn on the plasmid map.
    pub plasmid_show_features: bool,
    /// Whether restriction enzyme ticks are drawn on the plasmid map.
    pub plasmid_show_re: bool,
    /// Selected feature index in the plasmid feature list panel.
    pub plasmid_selected_feat: Option<usize>,
    /// Edit buffer: feature name.
    pub plasmid_edit_name: String,
    /// Edit buffer: feature start (1-based string for display).
    pub plasmid_edit_start: String,
    /// Edit buffer: feature end (1-based string).
    pub plasmid_edit_end: String,
    /// Edit buffer: feature color (buffered while editing, applied on Apply).
    pub plasmid_edit_color: helixview_core::color::Color,
    /// Edit buffer: hex string "#RRGGBB" for the color field.
    pub plasmid_edit_hex: String,
    /// Zoom factor for the plasmid canvas (1.0 = default, 0.2–5.0).
    pub plasmid_zoom: f32,
    /// Pan offset applied to the plasmid canvas center (pixels).
    pub plasmid_pan_x: f32,
    pub plasmid_pan_y: f32,
    /// Pre-computed RE sites for the active plasmid sequence.
    /// Keyed on (alignment ptr, seq_idx) so it's recomputed only when the sequence changes.
    pub plasmid_re_cache: Option<Arc<Vec<helixview_analysis::RestrictionSite>>>,
    /// Whether to show the covariation pairing-arcs strip above the entropy ruler.
    pub show_pairing_arcs: bool,
    /// Minimum covariation score (cov_score) for arcs to be displayed.
    pub pairing_arc_threshold: f64,
    /// Whether the hydrophobicity profile view is open.
    pub show_hydrophobicity: bool,
    /// Window size for Kyte-Doolittle hydrophobicity profile.
    pub hydro_window: usize,
    /// Whether the text alignment export view is open.
    pub show_text_export: bool,
    /// Options for text alignment export.
    pub text_export_opts: crate::text_export::TextExportOptions,
    /// Whether the command palette overlay is open.
    pub show_command_palette: bool,
    /// Current query string in the command palette.
    pub palette_query: String,
    /// Currently highlighted row in the command palette.
    pub palette_selected: usize,
    /// Background computation in progress (shown in status bar).
    pub computing: Option<String>,
    /// Pre-computed identity matrix (None = not yet computed or invalidated).
    pub identity_matrix_data: Option<Arc<Vec<Vec<f64>>>>,
    /// Whether the mutual information view is open.
    pub show_mutual_info: bool,
    /// Most recently computed MI result.
    pub mi_result: Option<Arc<helixview_analysis::MiResult>>,
    /// Whether MI computation is in progress.
    pub mi_running: bool,
    /// Minimum observations for MI computation.
    pub mi_min_obs: usize,
    /// Maximum top pairs to retain.
    pub mi_top_n: usize,
    /// Minimum distinct WC pair types for stem highlighting (1–6).
    pub mi_stem_min_wc_types: u8,
    /// Minimum WC fraction for stem highlighting (0–1).
    pub mi_stem_min_wc_frac: f64,
    /// Arc pointer of the alignment when MI was last run (for stale detection).
    pub mi_aln_ptr: Option<usize>,
    /// Active sub-tab in the MI view.
    pub mi_tab: crate::views::MiTab,
    /// User-configured accessory apps.
    pub accessories: Vec<AccessoryApp>,
    /// Whether the accessory manager view is open.
    pub show_accessories: bool,
    /// Index of the accessory being edited (None = new).
    pub editing_accessory: Option<usize>,
    /// Edit buffer: name field.
    pub acc_edit_name: String,
    /// Edit buffer: command field.
    pub acc_edit_cmd: String,
}

impl std::fmt::Debug for HelixViewApp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HelixViewApp")
            .field("status", &self.status)
            .field("seq_editor_idx", &self.seq_editor_idx)
            .finish_non_exhaustive()
    }
}

impl Default for HelixViewApp {
    fn default() -> Self {
        Self {
            document:        None,
            status:          String::new(),
            view:            ViewState::default(),
            history:         EditHistory::new(100),
            pairwise_result:       None,
            orf_results:           None,
            show_identity_matrix:  false,
            show_col_summary:      false,
            title_menu:            None,
            rename_input:          String::new(),
            show_restr_map:        false,
            show_dot_plot:         false,
            seq_editor_idx:        None,
            seq_editor_content:    None,
            color_table:           Arc::new(ColorTable::default()),
            show_color_editor:      false,
            color_editor_selected:  None,
            show_conservation:      false,
            conservation_threshold: 0.80,
            conservation_min_width: 4,
            six_frame_result:       None,
            protein_view_doc:       None,
            protein_view_frame:     0,
            blast:                  BlastState::Idle,
            dot_plot_seq_a:         0,
            dot_plot_seq_b:        1,
            dot_plot_window:       7,
            dot_plot_min_match:    5,
            show_composition:       false,
            show_oligo_tm:          false,
            tabs:                   vec![TabRecord::empty()],
            active_tab:             0,
            tree_data:              None,
            tree_selected:          BTreeSet::new(),
            trace_data:             None,
            trace_scroll:           0.0,
            trace_zoom:             1.0,
            trace_quant_ch_a:       b'A',
            trace_quant_ch_b:       b'C',
            prefs:                  crate::prefs::load(),
            show_prefs:             false,
            show_help:              false,
            show_shaded_export:     false,
            shaded_export_opts:     ShadedExportOptions::default(),
            show_taxonomy:          false,
            show_plasmid:           false,
            plasmid_seq_idx:        0,
            plasmid_show_features:  true,
            plasmid_show_re:        true,
            plasmid_selected_feat:  None,
            plasmid_edit_name:      String::new(),
            plasmid_edit_start:     String::new(),
            plasmid_edit_end:       String::new(),
            plasmid_edit_color:     helixview_core::color::Color::rgb(0.10, 0.60, 0.10),
            plasmid_edit_hex:       "#1a9919".to_string(),
            plasmid_zoom:           1.0,
            plasmid_pan_x:          0.0,
            plasmid_pan_y:          0.0,
            plasmid_re_cache:       None,
            computing:              None,
            identity_matrix_data:   None,
            show_pairing_arcs:      false,
            pairing_arc_threshold:  2.0,
            show_hydrophobicity:    false,
            hydro_window:           9,
            show_text_export:       false,
            text_export_opts:       crate::text_export::TextExportOptions::default(),
            show_command_palette:   false,
            palette_query:          String::new(),
            palette_selected:       0,
            show_mutual_info:       false,
            mi_result:              None,
            mi_running:             false,
            mi_min_obs:             2,
            mi_top_n:               100,
            mi_stem_min_wc_types:   2,
            mi_stem_min_wc_frac:    0.5,
            mi_aln_ptr:             None,
            mi_tab:                 crate::views::MiTab::TopPairs,
            accessories:            crate::accessories::load_accessories(),
            show_accessories:       false,
            editing_accessory:      None,
            acc_edit_name:          String::new(),
            acc_edit_cmd:           String::new(),
        }
    }
}

/// Which category of toolbar buttons is shown in the ribbon's second row.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum ToolbarTab {
    /// Sequence editing operations.
    #[default]
    Edit,
    /// Analysis and viewers.
    Analysis,
    /// External aligners and color / view options.
    View,
}

/// Whether the search bar is looking through residue patterns or sequence names.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum SearchMode {
    /// Search for a residue / IUPAC pattern across all sequences.
    #[default]
    Sequence,
    /// Search for a substring in sequence names / titles.
    Name,
}

/// How residue cells are colored in the grid.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum ColorScheme {
    /// Color every residue by its type (nucleotide or amino acid family).
    #[default]
    ResidueType,
    /// Color only residues that match the column consensus at ≥ 50% identity.
    Identity,
    /// No background coloring — plain white with dark text.
    Plain,
    /// Color by residue type, faded toward white for variable columns.
    Strength,
}

/// All per-document viewport state (scroll position, selection).
/// Reset when a new document is opened.
#[derive(Debug, Clone)]
pub struct ViewState {
    pub scroll_col:   usize,
    pub scroll_row:   usize,
    pub selected:     BTreeSet<usize>,
    pub color_scheme: ColorScheme,
    pub hover_col:    Option<usize>,
    pub hover_row:    Option<usize>,
    pub ctrl_held:     bool,
    pub shift_held:    bool,
    pub last_clicked:  Option<usize>,
    pub show_analysis: bool,
    /// Zoom level (1.0 = default, range 0.4–3.0).
    pub zoom:          f32,
    /// Whether the search bar is visible.
    pub show_search:   bool,
    /// Current text in the search bar.
    pub search_query:  String,
    /// Whether the search bar is in residue-pattern or name mode.
    pub search_mode:   SearchMode,
    /// Single-cell selection for residue editing (row, col).
    pub edit_cell:     Option<(usize, usize)>,
    /// Columns selected in the ruler (for bulk delete).
    pub selected_cols: BTreeSet<usize>,
    /// Whether to draw feature annotation stripes below each sequence row.
    pub show_features: bool,
    /// Residue/IUPAC pattern search results: list of (seq_idx, col) hits.
    pub pattern_matches: Vec<(usize, usize)>,
    /// Index into pattern_matches pointing to the currently highlighted hit.
    pub pattern_match_idx: usize,
    /// Whether the NCBI Entrez fetch bar is visible.
    pub show_fetch_bar: bool,
    /// Current text in the accession input field.
    pub fetch_input: String,
    /// Whether the annotation bar is visible (to annotate selected columns).
    pub show_annotate_bar: bool,
    /// Name field for the pending annotation.
    pub annotate_name: String,
    /// Feature type string for the pending annotation.
    pub annotate_type: String,
    /// When Some(frame), the protein view is active for this reading frame (0-based 0–5).
    pub protein_view_frame: Option<usize>,
    /// Whether the BLAST panel is open.
    pub show_blast_panel: bool,
    /// Program selected in the BLAST panel ("blastn"/"blastp"/"blastx").
    pub blast_program: String,
    /// Database selected in the BLAST panel.
    pub blast_database: String,
    /// Which category is active in the toolbar ribbon second row.
    pub toolbar_tab: ToolbarTab,
}

impl Default for ViewState {
    fn default() -> Self {
        Self {
            scroll_col:    0,
            scroll_row:    0,
            selected:      Default::default(),
            color_scheme:  ColorScheme::default(),
            hover_col:     None,
            hover_row:     None,
            ctrl_held:     false,
            shift_held:    false,
            last_clicked:  None,
            show_analysis: false,
            zoom:          1.0,
            show_search:   false,
            search_query:  String::new(),
            search_mode:   SearchMode::default(),
            edit_cell:     None,
            selected_cols:      Default::default(),
            show_features:      true,
            pattern_matches:    Vec::new(),
            pattern_match_idx:  0,
            show_fetch_bar:     false,
            fetch_input:        String::new(),
            show_annotate_bar:  false,
            annotate_name:      String::new(),
            annotate_type:      "misc_feature".to_string(),
            protein_view_frame: None,
            show_blast_panel:   false,
            blast_program:      "blastn".to_string(),
            blast_database:     "nt".to_string(),
            toolbar_tab:        ToolbarTab::default(),
        }
    }
}

// ── Messages ──────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub enum Message {
    /// User pressed "Open File…"
    OpenFileDialog,
    /// A file was chosen by the native open dialog.
    FileChosen(Option<std::path::PathBuf>),
    /// The file has been parsed (may contain an error string).
    FileLoaded(Result<SeqAlignment, String>),
    /// User pressed "Save As…"
    SaveFileDialog,
    /// A file path was chosen by the native save dialog.
    SaveChosen(Option<std::path::PathBuf>),
    /// User pressed "New Alignment"
    NewAlignment,
    /// Canvas grid reported a scroll delta (positive dy = scroll down).
    GridScrolled { dx: isize, dy: isize },
    /// User clicked a sequence row in the title panel.
    SequenceClicked { idx: usize, ctrl: bool, shift: bool },
    /// Keyboard modifier state changed.
    ModifiersChanged(Modifiers),
    /// Horizontal scrollbar slider moved.
    ScrollColSet(u32),
    /// User changed the color scheme via the toolbar.
    SetColorScheme(ColorScheme),
    /// Canvas grid reported current hover position.
    HoverPosition { col: Option<usize>, row: Option<usize> },
    /// Right-click on a grid column: insert gap if residue column, delete if all-gap column.
    GridRightClick { col: usize },
    /// User dragged a sequence row from `from` to `to`.
    SequenceMoved { from: usize, to: usize },
    /// Undo last edit (Ctrl+Z).
    Undo,
    /// Redo last undone edit (Ctrl+Y / Ctrl+Shift+Z).
    Redo,
    /// Scroll to the last column.
    ScrollToEnd,
    /// Select all sequences.
    SelectAll,
    /// Deselect all sequences.
    SelectNone,
    /// Delete gap columns near the hovered column (or all-gap columns in selection).
    DeleteGapColumns,
    /// Increase zoom level.
    ZoomIn,
    /// Decrease zoom level.
    ZoomOut,
    /// Reset zoom to 1.0.
    ZoomReset,
    /// User click-selected a column range in the ruler.
    ColumnsSelected { start: usize, end: usize },
    /// Clear column selection.
    ClearColumnSelection,
    /// Delete all selected columns (using history).
    DeleteSelectedColumns,
    /// Insert a gap at the cursor position in the active sequence row (Space).
    InsertGapAtCursor,
    /// Insert a gap at the cursor position in all OTHER sequence rows (Shift+Space).
    InsertGapInOthers,
    /// Delete the gap character immediately before the cursor in the active row (Backspace in edit mode).
    DeleteGapAtCursor,
    /// Delete the gap character immediately before the cursor in all OTHER rows (Shift+Backspace in edit mode).
    DeleteGapInOthers,
    /// Space key pressed; routed to gap insertion or ignored based on edit mode.
    SpacePressed { shift: bool },
    /// Backspace key pressed; routed to gap deletion or column deletion based on edit mode.
    BackspacePressed { shift: bool },
    /// Run ORF finder on the first selected/hovered sequence.
    RunOrfFinder,
    /// Close the ORF finder results panel.
    CloseOrfs,
    /// Run six-frame translation on the first selected/hovered sequence.
    RunSixFrame,
    /// Close the six-frame translation panel.
    CloseSixFrame,
    /// Toggle the analysis side panel on/off.
    ToggleAnalysis,
    /// Run NW pairwise alignment on the two selected sequences.
    RunPairwise,
    /// Close the pairwise result panel.
    ClosePairwise,
    /// Remove all all-gap columns from the alignment.
    MinimizeAlignment,
    /// Sort sequences alphabetically by name.
    SortByName,
    /// Reverse-complement the selected sequence rows (DNA only).
    ReverseComplementSelected,
    /// Toggle the sequence name search bar (Ctrl+F).
    ToggleSearch,
    /// User typed in the search bar.
    SearchInput(String),
    /// User pressed Escape to close the search bar.
    SearchClose,
    /// User clicked a residue cell (for single-cell editing).
    CellClicked { row: usize, col: usize },
    /// User typed a residue key while a cell is selected.
    ResidueTyped(char),
    /// Translate selected sequences (DNA→protein) in place.
    TranslateSelected,
    /// Append a consensus sequence row to the alignment.
    AddConsensus,
    /// Open the identity matrix view.
    ShowIdentityMatrix,
    /// Close the identity matrix view.
    CloseIdentityMatrix,
    /// Toggle feature stripe rendering below sequence rows.
    ToggleFeatures,
    /// Switch the search bar between Sequence and Name modes.
    SetSearchMode(SearchMode),
    /// Search for a residue pattern (IUPAC / amino acid) across all sequences.
    SearchPattern(String),
    /// Search for a substring across all sequence names.
    SearchNames(String),
    /// Jump to next/previous pattern match: +1 or -1.
    PatternMatchJump(isize),
    /// User right-clicked a sequence title label.
    TitleRightClick { idx: usize },
    /// Close the title context menu without action.
    TitleMenuClose,
    /// Rename the sequence at `idx` to a new name.
    RenameSequence { idx: usize, name: String },
    /// Toggle the lock state of the sequence at `idx`.
    ToggleLockSequence { idx: usize },
    /// Duplicate the sequence at `idx` (insert a copy right after).
    DuplicateSequence { idx: usize },
    /// Delete the sequence at `idx`.
    DeleteSequence { idx: usize },
    /// Update the rename text field value.
    RenameInput(String),
    /// Open the positional column summary view.
    ShowColSummary,
    /// Close the positional column summary view.
    CloseColSummary,
    /// Toggle the restriction map strip.
    ToggleRestrMap,
    /// Open the color table editor.
    OpenColorEditor,
    /// Close the color table editor.
    CloseColorEditor,
    /// Select a residue for RGB editing: (byte, is_nuc).
    SelectColorResidue(u8, bool),
    /// Adjust one channel of the currently selected residue color.
    /// `channel`: 0=R, 1=G, 2=B.
    SetResidueColor { residue: u8, is_nuc: bool, channel: u8, value: f32 },
    /// Reset the color table to built-in defaults.
    ResetColorTable,
    /// Open the conservation region search view.
    ShowConservation,
    /// Close the conservation region search view.
    CloseConservation,
    /// Update the conservation identity threshold.
    ConservationThreshold(f64),
    /// Update the conservation minimum window width.
    ConservationMinWidth(usize),
    /// Scroll the alignment to a specific column (used by search results).
    JumpToColumn(usize),
    /// Run external multiple-sequence aligner (all sequences).
    RunExternalAlign(ExternalAligner),
    /// Alignment result returned from external aligner.
    ExternalAlignDone(Result<SeqAlignment, String>),
    /// Open the raw-sequence editor for the sequence at `idx` (double-click title).
    EditSequenceRaw { idx: usize },
    /// The text_editor inside the sequence editor emitted an action.
    SeqEditorAction(text_editor::Action),
    /// Commit the edited raw sequence back to the alignment.
    SeqEditorCommit,
    /// Close the sequence editor without saving.
    SeqEditorClose,
    /// Open the dot-plot view.
    ShowDotPlot,
    /// Close the dot-plot view.
    CloseDotPlot,
    /// Change the X-axis sequence index for the dot plot.
    DotPlotSeqA(usize),
    /// Change the Y-axis sequence index for the dot plot.
    DotPlotSeqB(usize),
    /// Change the dot-plot window size.
    DotPlotWindow(usize),
    /// Change the dot-plot minimum-match threshold.
    DotPlotMinMatch(usize),

    // ── NCBI Entrez fetch ────────────────────────────────────────────────────
    /// Show/hide the accession fetch bar.
    ToggleFetchBar,
    /// User is typing in the accession field.
    FetchInput(String),
    /// Trigger an NCBI fetch for the current accession input.
    FetchAccession,
    /// Result returned from the background fetch task.
    EntrezResult(Result<SeqAlignment, String>),

    // ── Export ───────────────────────────────────────────────────────────────
    /// Open a save dialog and export selected sequences (or all) as FASTA.
    ExportFastaDialog,
    /// Export gap-stripped sequences as FASTA.
    ExportRawFastaDialog,
    /// A save path was chosen for FASTA export; `raw` strips gaps.
    ExportFastaChosen { path: Option<std::path::PathBuf>, raw: bool },
    /// Export the identity matrix as CSV.
    ExportIdentityCsv,
    /// A save path was chosen for the CSV export.
    ExportCsvChosen(Option<std::path::PathBuf>),

    // ── Feature annotation ───────────────────────────────────────────────────
    /// Show/hide the inline annotation bar.
    ToggleAnnotateBar,
    /// User typed in the feature name field.
    AnnotateName(String),
    /// User changed the feature type string.
    AnnotateType(String),
    /// Apply the annotation to selected sequences over selected columns.
    ApplyAnnotation,

    // ── Protein view / translation toggle ────────────────────────────────────
    /// Enter or refresh protein view for reading frame `frame` (0-based 0–5).
    SetProteinViewFrame(usize),
    /// Return to DNA view (clear protein_view_doc).
    ExitProteinView,

    // ── BLAST ─────────────────────────────────────────────────────────────────
    /// Open/close the BLAST panel.
    ToggleBlastPanel,
    /// User changed the BLAST program.
    BlastProgramChanged(String),
    /// User changed the BLAST database.
    BlastDatabaseChanged(String),
    /// Submit the selected/hovered sequence to NCBI BLAST.
    SubmitBlast,
    /// RID received from NCBI (submission succeeded).
    BlastSubmitted(Result<String, String>),
    /// User clicked "Check" to poll for results.
    BlastCheck,
    /// Poll result: Ok(true) = ready, Ok(false) = still waiting.
    BlastPollResult(Result<bool, String>),
    /// Full results fetched from NCBI.
    BlastResultsFetched(Result<Vec<BlastHit>, String>),
    /// Import selected BLAST hits by accession via Entrez fetch.
    BlastImportHit(String),
    /// Close the BLAST results view.
    CloseBlast,

    // ── Toolbar ribbon ────────────────────────────────────────────────────────
    /// Switch the ribbon second row to show a different button category.
    SetToolbarTab(ToolbarTab),

    // ── Composition & Tm ─────────────────────────────────────────────────────
    ShowComposition,
    CloseComposition,
    ShowOligoTm,
    CloseOligoTm,

    // ── Hydrophobicity profile ────────────────────────────────────────────────
    ShowHydrophobicity,
    CloseHydrophobicity,
    HydroWindow(usize),

    // ── PDF export ────────────────────────────────────────────────────────────
    ExportPdfDialog,
    PdfExportPath(Option<std::path::PathBuf>),
    PdfExportDone(Result<(), String>),

    // ── Text alignment export ─────────────────────────────────────────────────
    ShowTextExport,
    CloseTextExport,
    TextExportResPerRow(usize),
    TextExportTitleChars(usize),
    TextExportShowRuler(bool),
    TextExportNumberLines(bool),
    TextExportSave,
    TextExportPath(Option<std::path::PathBuf>),
    TextExportDone(Result<(), String>),

    // ── Async identity matrix ─────────────────────────────────────────────────
    ComputeIdentityMatrix,
    IdentityMatrixDone(Arc<Vec<Vec<f64>>>),

    // ── Async pairwise ────────────────────────────────────────────────────────
    PairwiseDone(Result<helixview_analysis::PairwiseResult, String>),

    // ── Local BLAST DB ────────────────────────────────────────────────────────
    CreateBlastDbDialog,
    CreateBlastDbFasta(Option<std::path::PathBuf>),
    CreateBlastDbDone(Result<String, String>),

    // ── Accessory apps ────────────────────────────────────────────────────────
    ShowAccessories,
    CloseAccessories,
    RunAccessory(usize),
    AccessoryDone(Result<String, String>),
    EditAccessory(usize),
    SaveAccessory,
    DeleteAccessory(usize),
    NewAccessory,
    AccessoryFieldName(String),
    AccessoryFieldCmd(String),

    // ── Mutual information ────────────────────────────────────────────────────
    ShowMutualInfo,
    CloseMutualInfo,
    RunMutualInfo,
    MiDone(Result<Arc<helixview_analysis::MiResult>, String>),
    MiMinObs(usize),
    MiTopN(usize),
    MiStemMinWcTypes(u8),
    MiStemMinWcFrac(f64),
    MiSetTab(crate::views::MiTab),
    /// Toggle the covariation pairing-arcs strip in the alignment view.
    TogglePairingArcs,
    /// Set the minimum covariation score threshold for displayed arcs.
    PairingArcThreshold(f64),

    // ── Command palette ───────────────────────────────────────────────────────
    OpenCommandPalette,
    CloseCommandPalette,
    PaletteQuery(String),
    PaletteSelect(usize),
    PaletteConfirm,
    PaletteUp,
    PaletteDown,

    // ── Help / shortcuts ──────────────────────────────────────────────────────
    ShowHelp,
    CloseHelp,

    // ── Project file ──────────────────────────────────────────────────────────
    SaveProjectDialog,
    SaveProjectChosen(Option<std::path::PathBuf>),
    SaveProjectDone(Result<(), String>),
    OpenProjectDialog,
    OpenProjectChosen(Option<std::path::PathBuf>),
    OpenProjectDone(Result<(), String>),

    // ── Preferences ───────────────────────────────────────────────────────────
    OpenPrefs,
    ClosePrefs,
    SavePrefs,
    PrefsUndoDepth(usize),
    PrefsDefaultZoom(f32),
    PrefsFontSize(f32),
    PrefsShowFeatures(bool),
    PrefsShowRestrMap(bool),
    PrefsIdentityThreshold(f32),
    PrefsRecentFilesMax(usize),
    PrefsClearRecent,
    OpenRecentFile(std::path::PathBuf),

    // ── Shaded graphic export ─────────────────────────────────────────────────
    /// Open the shaded export options dialog.
    ShowShadedExportDialog,
    /// Close the options dialog without exporting.
    CloseShadedExportDialog,
    /// Change the shading scheme in the export dialog.
    ShadedExportSetScheme(ColorScheme),
    /// Change the identity threshold slider value.
    ShadedExportThreshold(f32),
    /// Toggle ruler row.
    ShadedExportToggleRuler(bool),
    /// Toggle consensus row.
    ShadedExportToggleConsensus(bool),
    /// Open file-save dialog and run the export.
    ShadedExportRun,
    /// Save path chosen.
    ShadedExportChosen(Option<std::path::PathBuf>),
    /// Background write completed.
    ShadedExportDone(Result<(), String>),

    // ── Taxonomy table ────────────────────────────────────────────────────────
    ShowTaxonomy,
    CloseTaxonomy,
    /// Reorder alignment sequences by taxonomy rank level.
    /// `level`: 0-based rank index into the lineage; `usize::MAX` = sort by organism name.
    SortByTaxonomy(usize),

    // ── SVG export ────────────────────────────────────────────────────────────
    /// Open save dialog for SVG alignment export.
    ExportSvgDialog,
    /// Save path chosen for SVG export.
    ExportSvgChosen(Option<std::path::PathBuf>),
    /// Background SVG write completed.
    ExportSvgDone(Result<(), String>),

    // ── Tabs ──────────────────────────────────────────────────────────────────
    /// Create a new empty tab and switch to it.
    NewTab,
    /// Switch to the tab at this index.
    SwitchTab(usize),
    /// Close the tab at this index.
    CloseTab(usize),

    // ── Phylogenetic tree viewer ──────────────────────────────────────────────
    /// Build a NJ tree from the current alignment (no file picker).
    BuildTreeFromAlignment,
    /// Open a file-picker dialog for .nwk / .nex / .nexus / .tree files.
    OpenTreeDialog,
    /// A tree file path was chosen (None = cancelled).
    TreeFileChosen(Option<std::path::PathBuf>),
    /// Tree file was parsed (may be an error).
    TreeLoaded(Result<Vec<helixview_core::PhyloTree>, String>),
    /// Close the tree view and return to the alignment.
    CloseTree,
    /// User clicked a leaf node in the tree canvas.
    TreeNodeClicked(String),
    /// Clear the leaf selection in the tree viewer.
    TreeClearSelection,
    /// Reorder the alignment rows to match the current tree leaf order.
    TreeReorderAlignment,

    // ── Plasmid map viewer ────────────────────────────────────────────────────
    /// Open the circular plasmid map view for the first (or selected) sequence.
    ShowPlasmid,
    /// Close the plasmid map view.
    ClosePlasmid,
    /// Navigate to the previous sequence in the plasmid view.
    PlasmidPrevSeq,
    /// Navigate to the next sequence in the plasmid view.
    PlasmidNextSeq,
    /// Toggle feature arc visibility in the plasmid view.
    PlasmidToggleFeatures,
    /// Toggle restriction enzyme tick visibility in the plasmid view.
    PlasmidToggleRe,
    /// Background computation of plasmid RE sites finished.
    PlasmidReSitesDone(Vec<helixview_analysis::RestrictionSite>),
    /// Open a save dialog for plasmid SVG export.
    ExportPlasmidSvg,
    /// Save path chosen for plasmid SVG export.
    ExportPlasmidSvgChosen(Option<std::path::PathBuf>),
    /// Background plasmid SVG write completed.
    ExportPlasmidSvgDone(Result<(), String>),
    /// Select (or deselect) a feature in the plasmid feature list panel.
    PlasmidSelectFeature(Option<usize>),
    /// Edit buffer changes for plasmid feature editing.
    PlasmidEditName(String),
    PlasmidEditStart(String),
    PlasmidEditEnd(String),
    PlasmidEditColor(helixview_core::color::Color),
    PlasmidEditHex(String),
    /// Commit edits to the selected feature.
    PlasmidApplyEdit,
    /// Delete the feature at the given index.
    PlasmidDeleteFeature(usize),
    /// Add a blank new feature to the current sequence.
    PlasmidAddFeature,
    /// Continuous zoom delta from mouse wheel (multiplicative).
    PlasmidZoomDelta(f32),
    /// Step zoom in.
    PlasmidZoomIn,
    /// Step zoom out.
    PlasmidZoomOut,
    /// Reset zoom to 1.0 and pan to (0, 0).
    PlasmidZoomReset,
    /// Drag-to-pan delta from canvas mouse drag.
    PlasmidPanDelta(f32, f32),

    // ── ABI trace viewer ──────────────────────────────────────────────────────
    /// Open a file-picker dialog for .ab1 / .abi trace files.
    OpenTraceDialog,
    /// A trace file path was chosen (None = cancelled).
    TraceFileChosen(Option<std::path::PathBuf>),
    /// Trace file was parsed and returned (may be an error).
    TraceLoaded(Result<helixview_formats::abi::AbiTrace, String>),
    /// Close the trace view and return to the alignment.
    CloseTrace,
    /// Set the horizontal scroll fraction (0.0–1.0) of the chromatogram.
    TraceScroll(f32),
    /// Set the zoom level of the chromatogram (1.0–16.0).
    TraceZoom(f32),
    /// Set quantification channel A to this base letter.
    TraceQuantChA(u8),
    /// Set quantification channel B to this base letter.
    TraceQuantChB(u8),

    // ── Application lifecycle ─────────────────────────────────────────────────
    /// Window close was requested; save session and exit.
    AppCloseRequested(iced::window::Id),
}

// ── Tab helpers ───────────────────────────────────────────────────────────────

impl HelixViewApp {
    /// Display name for a tab (document name, or "New Tab").
    pub fn tab_name(rec: &TabRecord) -> String {
        rec.document.as_ref()
            .map(|d| d.name.clone())
            .filter(|n| !n.is_empty())
            .unwrap_or_else(|| "New Tab".to_string())
    }

    /// Snapshot the current live state into `tabs[active_tab]` and switch to `new_idx`.
    /// Clears all transient analysis state. Does nothing if `new_idx == active_tab`.
    pub fn switch_to(&mut self, new_idx: usize) {
        if new_idx == self.active_tab || new_idx >= self.tabs.len() { return; }

        // Save live state into the current slot.
        self.tabs[self.active_tab] = TabRecord {
            document:           self.document.clone(),
            view:               self.view.clone(),
            history:            std::mem::replace(&mut self.history, EditHistory::new(100)),
            show_restr_map:     self.show_restr_map,
            protein_view_doc:   self.protein_view_doc.clone(),
            protein_view_frame: self.protein_view_frame,
        };

        // Restore from new slot.
        let rec = std::mem::replace(&mut self.tabs[new_idx], TabRecord::empty());
        self.document           = rec.document;
        self.view               = rec.view;
        self.history            = rec.history;
        self.show_restr_map     = rec.show_restr_map;
        self.protein_view_doc   = rec.protein_view_doc;
        self.protein_view_frame = rec.protein_view_frame;
        self.active_tab         = new_idx;

        // Clear transient analysis state.
        self.pairwise_result    = None;
        self.orf_results        = None;
        self.six_frame_result   = None;
        self.blast              = BlastState::Idle;
        self.trace_data         = None;
        self.trace_scroll       = 0.0;
        self.trace_zoom         = 1.0;
        self.tree_data          = None;
        self.tree_selected.clear();
        self.show_identity_matrix    = false;
        self.identity_matrix_data    = None;
        self.show_col_summary   = false;
        self.show_dot_plot      = false;
        self.show_conservation  = false;
        self.seq_editor_idx     = None;
        self.seq_editor_content = None;
        self.title_menu         = None;
    }

    /// Open a document in a tab: use current tab if empty, else create a new tab.
    /// Applies user preferences (zoom, features, restr map) to the new view.
    pub fn open_in_tab(&mut self, aln: SeqAlignment) {
        let current_is_empty = self.document.is_none();
        if current_is_empty {
            self.document = Some(Arc::new(aln));
            let mut v = ViewState::default();
            v.zoom = self.prefs.default_zoom;
            v.show_features = self.prefs.show_features;
            self.view = v;
            self.show_restr_map = self.prefs.show_restr_map;
            self.history = EditHistory::new(self.prefs.undo_depth);
        } else {
            // Snapshot current and open new tab.
            self.tabs[self.active_tab] = TabRecord {
                document:           self.document.clone(),
                view:               self.view.clone(),
                history:            std::mem::replace(&mut self.history, EditHistory::new(100)),
                show_restr_map:     self.show_restr_map,
                protein_view_doc:   self.protein_view_doc.clone(),
                protein_view_frame: self.protein_view_frame,
            };
            let new_idx = self.tabs.len();
            self.tabs.push(TabRecord::empty());
            self.active_tab = new_idx;
            self.document = Some(Arc::new(aln));
            let mut v = ViewState::default();
            v.zoom = self.prefs.default_zoom;
            v.show_features = self.prefs.show_features;
            self.view = v;
            self.history = EditHistory::new(self.prefs.undo_depth);
            self.show_restr_map = self.prefs.show_restr_map;
            self.protein_view_doc = None;
            self.protein_view_frame = 0;
            self.pairwise_result  = None;
            self.orf_results      = None;
            self.six_frame_result = None;
            self.blast            = BlastState::Idle;
            self.trace_data       = None;
            self.tree_data        = None;
            self.tree_selected.clear();
        }
    }
}

// ── Update ────────────────────────────────────────────────────────────────────

pub fn update(app: &mut HelixViewApp, message: Message) -> Task<Message> {
    match message {
        Message::OpenFileDialog => {
            Task::perform(pick_file(), Message::FileChosen)
        }

        Message::FileChosen(None) => Task::none(),

        Message::FileChosen(Some(path)) => {
            crate::prefs::push_recent(&mut app.prefs, path.clone());
            let _ = crate::prefs::save(&app.prefs);
            app.status = format!("Loading {}…", path.display());
            Task::perform(load_file(path), Message::FileLoaded)
        }

        Message::FileLoaded(Ok(aln)) => {
            app.status = format!(
                "Loaded {} sequences, {} columns",
                aln.seq_count(), aln.col_count()
            );
            app.open_in_tab(aln);
            app.plasmid_re_cache = None;
            auto_save_session(app);
            Task::none()
        }

        Message::FileLoaded(Err(e)) => {
            app.status = format!("Error: {e}");
            Task::none()
        }

        Message::SaveFileDialog => {
            if app.document.is_some() {
                Task::perform(save_file_dialog(), Message::SaveChosen)
            } else {
                Task::none()
            }
        }

        Message::SaveChosen(None) => Task::none(),

        Message::SaveChosen(Some(path)) => {
            if let Some(aln) = &app.document {
                match helixview_formats::save(aln, &path) {
                    Ok(()) => app.status = format!("Saved to {}", path.display()),
                    Err(e) => app.status = format!("Save error: {e}"),
                }
            }
            Task::none()
        }

        Message::NewAlignment => {
            app.open_in_tab(SeqAlignment::new("Untitled"));
            app.status = "New alignment created".to_string();
            Task::none()
        }

        Message::GridScrolled { dx, dy } => {
            if let Some(aln) = &app.document {
                let ncols = aln.col_count();
                let nrows = aln.seq_count();

                // When a cell is selected for editing, arrow keys move the cursor.
                if let Some((row, col)) = app.view.edit_cell {
                    let new_col = ((col as isize + dx).max(0) as usize).min(ncols.saturating_sub(1));
                    let new_row = ((row as isize + dy).max(0) as usize).min(nrows.saturating_sub(1));
                    app.view.edit_cell = Some((new_row, new_col));
                    // Auto-scroll to keep the edit cell visible (approx 100 cols / 30 rows viewport).
                    if new_col < app.view.scroll_col {
                        app.view.scroll_col = new_col;
                    } else if new_col >= app.view.scroll_col + 100 {
                        app.view.scroll_col = new_col.saturating_sub(99);
                    }
                    if new_row < app.view.scroll_row {
                        app.view.scroll_row = new_row;
                    } else if new_row >= app.view.scroll_row + 30 {
                        app.view.scroll_row = new_row.saturating_sub(29);
                    }
                    return Task::none();
                }

                app.view.scroll_col = clamp_add(app.view.scroll_col, dx, ncols);
                app.view.scroll_row = clamp_add(app.view.scroll_row, dy, nrows);
            }
            Task::none()
        }

        Message::SequenceClicked { idx, ctrl, shift } => {
            if shift {
                // Range select from last_clicked to idx (inclusive)
                let from = app.view.last_clicked.unwrap_or(idx);
                let (lo, hi) = if from <= idx { (from, idx) } else { (idx, from) };
                if !ctrl {
                    app.view.selected.clear();
                }
                for i in lo..=hi {
                    app.view.selected.insert(i);
                }
            } else if ctrl {
                // Toggle
                if !app.view.selected.remove(&idx) {
                    app.view.selected.insert(idx);
                }
            } else {
                // Single select
                app.view.selected.clear();
                app.view.selected.insert(idx);
            }
            app.view.last_clicked = Some(idx);
            Task::none()
        }

        Message::ModifiersChanged(mods) => {
            app.view.ctrl_held  = mods.control();
            app.view.shift_held = mods.shift();
            Task::none()
        }

        Message::ScrollColSet(col) => {
            app.view.scroll_col = col as usize;
            Task::none()
        }

        Message::SetColorScheme(scheme) => {
            app.view.color_scheme = scheme;
            Task::none()
        }

        Message::HoverPosition { col, row } => {
            app.view.hover_col = col;
            app.view.hover_row = row;
            Task::none()
        }

        Message::GridRightClick { col } => {
            if let Some(arc) = &mut app.document {
                let mut aln = (**arc).clone();
                if aln.is_gap_column(col) {
                    let cmd = Box::new(DelColumn::new(col));
                    app.history.execute(cmd, &mut aln);
                    app.status = format!(
                        "Deleted gap column {}  (Ctrl+Z to undo)",
                        col + 1,
                    );
                } else {
                    let cmd = Box::new(InsertGapColumn { col });
                    app.history.execute(cmd, &mut aln);
                    app.status = format!(
                        "Inserted gap at column {}  (Ctrl+Z to undo)",
                        col + 1,
                    );
                }
                *arc = Arc::new(aln);
            }
            Task::none()
        }

        Message::ScrollToEnd => {
            if let Some(aln) = &app.document {
                app.view.scroll_col = aln.col_count().saturating_sub(1);
            }
            Task::none()
        }

        Message::SelectAll => {
            if let Some(aln) = &app.document {
                app.view.selected = (0..aln.seq_count()).collect();
            }
            Task::none()
        }

        Message::SelectNone => {
            app.view.selected.clear();
            app.view.edit_cell = None;
            Task::none()
        }

        Message::DeleteGapColumns => {
            if let Some(arc) = &mut app.document {
                let mut aln = (**arc).clone();
                let before = aln.col_count();
                // Delete all all-gap columns (using minimize is the simplest approach).
                aln.minimize();
                let after   = aln.col_count();
                let removed = before.saturating_sub(after);
                if removed > 0 {
                    app.history = EditHistory::new(100);
                    *arc = Arc::new(aln);
                    app.status = format!("Removed {removed} gap column(s), {after} remaining");
                } else {
                    app.status = "No all-gap columns to remove".to_string();
                }
            }
            Task::none()
        }

        Message::ZoomIn => {
            app.view.zoom = (app.view.zoom + 0.25).min(3.0);
            Task::none()
        }
        Message::ZoomOut => {
            app.view.zoom = (app.view.zoom - 0.25).max(0.4);
            Task::none()
        }
        Message::ZoomReset => {
            app.view.zoom = 1.0;
            Task::none()
        }

        Message::ColumnsSelected { start, end } => {
            let (lo, hi) = if start <= end { (start, end) } else { (end, start) };
            app.view.selected_cols = (lo..=hi).collect();
            app.status = format!(
                "Selected {} column(s) ({}–{}). Backspace to delete.",
                hi - lo + 1, lo + 1, hi + 1
            );
            Task::none()
        }

        Message::DeleteSelectedColumns => {
            if let Some(arc) = &mut app.document {
                let cols: Vec<usize> = {
                    let mut v: Vec<usize> = app.view.selected_cols.iter().copied().collect();
                    v.sort_unstable_by(|a, b| b.cmp(a)); // delete from right to left
                    v
                };
                if !cols.is_empty() {
                    let mut aln = (**arc).clone();
                    let count   = cols.len();
                    for col in &cols {
                        let cmd = Box::new(DelColumn::new(*col));
                        app.history.execute(cmd, &mut aln);
                    }
                    *arc = Arc::new(aln);
                    app.view.selected_cols.clear();
                    app.status = format!("Deleted {count} column(s)  (Ctrl+Z to undo)");
                }
            }
            Task::none()
        }

        Message::InsertGapAtCursor => {
            if let (Some((row, col)), Some(arc)) = (app.view.edit_cell, &mut app.document) {
                let mut aln = (**arc).clone();
                let cmd = Box::new(InsertGapInSeq { seq_idx: row, col, count: 1 });
                app.history.execute(cmd, &mut aln);
                *arc = Arc::new(aln);
                // Advance cursor one column to stay after the inserted gap.
                app.view.edit_cell = Some((row, col + 1));
                app.status = "Inserted gap  (Ctrl+Z to undo)".to_string();
            }
            Task::none()
        }

        Message::InsertGapInOthers => {
            if let (Some((row, col)), Some(arc)) = (app.view.edit_cell, &mut app.document) {
                let mut aln = (**arc).clone();
                let n = aln.seq_count();
                for other in 0..n {
                    if other == row { continue; }
                    let cmd = Box::new(InsertGapInSeq { seq_idx: other, col, count: 1 });
                    app.history.execute(cmd, &mut aln);
                }
                *arc = Arc::new(aln);
                app.status = "Inserted gap in other rows  (Ctrl+Z to undo)".to_string();
            }
            Task::none()
        }

        Message::DeleteGapAtCursor => {
            if let (Some((row, col)), Some(arc)) = (app.view.edit_cell, &mut app.document) {
                // Backspace: delete gap immediately *before* cursor (col - 1).
                if col == 0 { return Task::none(); }
                let del_col = col - 1;
                let is_gap = arc.sequences.get(row)
                    .and_then(|s| s.residues.get(del_col))
                    .map(|&b| helixview_core::sequence::is_gap(b))
                    .unwrap_or(false);
                if is_gap {
                    let mut aln = (**arc).clone();
                    let cmd = Box::new(DeleteGapInSeq::new(row, del_col, 1));
                    app.history.execute(cmd, &mut aln);
                    *arc = Arc::new(aln);
                    app.view.edit_cell = Some((row, del_col));
                    app.status = "Deleted gap  (Ctrl+Z to undo)".to_string();
                }
            }
            Task::none()
        }

        Message::DeleteGapInOthers => {
            if let (Some((row, col)), Some(arc)) = (app.view.edit_cell, &mut app.document) {
                if col == 0 { return Task::none(); }
                let del_col = col - 1;
                let mut aln = (**arc).clone();
                let n = aln.seq_count();
                let mut deleted = 0usize;
                for other in 0..n {
                    if other == row { continue; }
                    let is_gap = aln.sequences.get(other)
                        .and_then(|s| s.residues.get(del_col))
                        .map(|&b| helixview_core::sequence::is_gap(b))
                        .unwrap_or(false);
                    if is_gap {
                        let cmd = Box::new(DeleteGapInSeq::new(other, del_col, 1));
                        app.history.execute(cmd, &mut aln);
                        deleted += 1;
                    }
                }
                if deleted > 0 {
                    *arc = Arc::new(aln);
                    app.status = format!("Deleted gap in {deleted} other row(s)  (Ctrl+Z to undo)");
                }
            }
            Task::none()
        }

        Message::SpacePressed { shift } => {
            if app.view.edit_cell.is_some() {
                let msg = if shift { Message::InsertGapInOthers } else { Message::InsertGapAtCursor };
                update(app, msg)
            } else {
                Task::none()
            }
        }

        Message::BackspacePressed { shift } => {
            if app.view.edit_cell.is_some() {
                let msg = if shift { Message::DeleteGapInOthers } else { Message::DeleteGapAtCursor };
                update(app, msg)
            } else {
                update(app, Message::DeleteSelectedColumns)
            }
        }

        Message::RunOrfFinder => {
            if let Some(aln) = &app.document {
                let idx = app.view.selected.iter().next().copied()
                    .or(app.view.hover_row)
                    .unwrap_or(0);
                if idx < aln.seq_count() {
                    let seq  = &aln.sequences[idx];
                    let orfs = helixview_analysis::find_orfs(&seq.residues, 25);
                    let name = seq.name.clone();
                    app.status = format!("Found {} ORF(s) in {}", orfs.len(), name);
                    app.orf_results = Some((name, orfs));
                }
            }
            Task::none()
        }

        Message::CloseOrfs => {
            app.orf_results = None;
            Task::none()
        }

        Message::RunSixFrame => {
            if let Some(aln) = &app.document {
                let idx = app.view.selected.iter().next().copied()
                    .or(app.view.hover_row)
                    .unwrap_or(0);
                if idx < aln.seq_count() {
                    let seq  = &aln.sequences[idx];
                    let name = seq.name.clone();
                    let nt_len = seq.residues.iter()
                        .filter(|&&b| !helixview_core::sequence::is_gap(b))
                        .count();
                    let frames = helixview_analysis::six_frame_all(&seq.residues);
                    app.status = format!("6-frame translation for {name}  ({nt_len} nt)");
                    app.six_frame_result = Some((name, nt_len, frames));
                }
            }
            Task::none()
        }

        Message::CloseSixFrame => {
            app.six_frame_result = None;
            Task::none()
        }

        Message::ToggleAnalysis => {
            app.view.show_analysis = !app.view.show_analysis;
            Task::none()
        }

        Message::RunPairwise => {
            if let Some(aln) = &app.document {
                let selected: Vec<usize> = app.view.selected.iter().copied().collect();
                if selected.len() == 2 {
                    let seq_a = aln.sequences[selected[0]].residues
                        .iter().filter(|&&b| !matches!(b, b'-' | b'~' | b'.'))
                        .copied().collect::<Vec<u8>>();
                    let seq_b = aln.sequences[selected[1]].residues
                        .iter().filter(|&&b| !matches!(b, b'-' | b'~' | b'.'))
                        .copied().collect::<Vec<u8>>();
                    app.computing = Some("Running pairwise alignment…".to_string());
                    Task::perform(
                        async move {
                            tokio::task::spawn_blocking(move || {
                                helixview_analysis::needleman_wunsch(
                                    &seq_a, &seq_b,
                                    helixview_analysis::DEFAULT_MATCH,
                                    helixview_analysis::DEFAULT_MISMATCH,
                                    helixview_analysis::DEFAULT_GAP,
                                )
                            })
                            .await
                            .map_err(|e| e.to_string())
                        },
                        Message::PairwiseDone,
                    )
                } else {
                    app.status = "Select exactly 2 sequences to run pairwise alignment".to_string();
                    Task::none()
                }
            } else {
                Task::none()
            }
        }

        Message::ClosePairwise => {
            app.pairwise_result = None;
            Task::none()
        }

        Message::MinimizeAlignment => {
            if let Some(arc) = &mut app.document {
                let mut aln = (**arc).clone();
                let before = aln.col_count();
                aln.minimize();
                let after  = aln.col_count();
                let removed = before.saturating_sub(after);
                *arc = Arc::new(aln);
                // Clear history since minimize is a non-trivial structural change
                app.history = EditHistory::new(100);
                app.status = format!("Minimized: removed {removed} all-gap columns ({after} remaining)");
            }
            Task::none()
        }

        Message::SortByName => {
            if let Some(arc) = &mut app.document {
                let mut aln = (**arc).clone();
                aln.sequences.sort_by(|a, b| a.name.cmp(&b.name));
                *arc = Arc::new(aln);
                app.view.selected.clear();
                app.history = EditHistory::new(100);
                app.status = "Sorted sequences by name".to_string();
            }
            Task::none()
        }

        Message::ReverseComplementSelected => {
            if let Some(arc) = &mut app.document {
                if app.view.selected.is_empty() {
                    app.status = "Select sequence rows first, then use Rev-Comp".to_string();
                } else {
                    let mut aln = (**arc).clone();
                    let count = app.view.selected.len();
                    for &idx in &app.view.selected {
                        if idx < aln.sequences.len() {
                            aln.sequences[idx].reverse_complement();
                        }
                    }
                    *arc = Arc::new(aln);
                    app.status = format!("Reverse-complemented {count} sequence(s)");
                }
            }
            Task::none()
        }

        Message::TranslateSelected => {
            if let Some(arc) = &mut app.document {
                if app.view.selected.is_empty() {
                    app.status = "Select sequence rows first, then use Translate".to_string();
                } else {
                    let mut aln = (**arc).clone();
                    let count = app.view.selected.len();
                    for &idx in &app.view.selected {
                        if idx < aln.sequences.len() {
                            // Strip gaps, translate frame 0, replace residues.
                            let raw: Vec<u8> = aln.sequences[idx].residues.iter()
                                .filter(|&&b| !matches!(b, b'-' | b'~' | b'.'))
                                .copied().collect();
                            let protein = helixview_analysis::translate_frame(&raw, 0);
                            aln.sequences[idx].residues = protein;
                        }
                    }
                    *arc = Arc::new(aln);
                    app.history = EditHistory::new(100);
                    app.status = format!("Translated {count} sequence(s) to protein");
                }
            }
            Task::none()
        }

        Message::ToggleSearch => {
            app.view.show_search = !app.view.show_search;
            if !app.view.show_search {
                app.view.search_query.clear();
            }
            Task::none()
        }

        Message::SetSearchMode(mode) => {
            app.view.search_mode = mode;
            app.view.pattern_matches.clear();
            app.view.pattern_match_idx = 0;
            // Re-run the search in the new mode if there is already a query.
            let q = app.view.search_query.clone();
            if !q.is_empty() {
                return update(app, Message::SearchInput(q));
            }
            Task::none()
        }

        Message::SearchInput(q) => {
            app.view.search_query = q.clone();
            app.view.pattern_matches.clear();
            app.view.pattern_match_idx = 0;
            match app.view.search_mode {
                SearchMode::Name => {
                    return update(app, Message::SearchNames(q));
                }
                SearchMode::Sequence => {
                    if q.len() >= 2 {
                        return update(app, Message::SearchPattern(q));
                    }
                }
            }
            Task::none()
        }

        Message::SearchClose => {
            // Escape priority: command palette first, then cell-edit, then column selection, then search bar.
            if app.show_command_palette {
                app.show_command_palette = false;
                return Task::none();
            }
            if app.view.edit_cell.take().is_some() {
                // done
            } else if !app.view.selected_cols.is_empty() {
                return update(app, Message::ClearColumnSelection);
            } else {
                app.view.show_search  = false;
                app.view.search_query.clear();
                app.view.pattern_matches.clear();
                app.view.pattern_match_idx = 0;
            }
            Task::none()
        }

        Message::ClearColumnSelection => {
            app.view.selected_cols.clear();
            Task::none()
        }

        Message::CellClicked { row, col } => {
            app.view.edit_cell = Some((row, col));
            Task::none()
        }

        Message::ResidueTyped(ch) => {
            if let (Some((row, col)), Some(arc)) =
                (app.view.edit_cell, &mut app.document)
            {
                let b = (ch as u8).to_ascii_uppercase();
                // Only accept valid residue chars and gap
                if b.is_ascii_alphabetic() || b == b'-' || b == b'.' {
                    let mut aln = (**arc).clone();
                    if row < aln.seq_count() {
                        let seq = &aln.sequences[row];
                        if col < seq.residues.len() {
                            let old_byte = seq.residues[col];
                            if old_byte != b {
                                let cmd = Box::new(SetResidues {
                                    seq_idx:   row,
                                    position:  col,
                                    old_bytes: vec![old_byte],
                                    new_bytes: vec![b],
                                });
                                app.history.execute(cmd, &mut aln);
                                *arc = Arc::new(aln);
                            }
                        }
                        // Advance to next column automatically.
                        let new_col = col + 1;
                        if let Some(arc2) = &app.document {
                            if new_col < arc2.col_count() {
                                app.view.edit_cell = Some((row, new_col));
                            }
                        }
                    }
                }
            }
            Task::none()
        }

        Message::ShowIdentityMatrix => {
            app.show_identity_matrix = true;
            // Always recompute when opening.
            return update(app, Message::ComputeIdentityMatrix);
        }

        Message::CloseIdentityMatrix => {
            app.show_identity_matrix = false;
            Task::none()
        }

        Message::ShowColSummary => {
            app.show_col_summary = true;
            Task::none()
        }

        Message::CloseColSummary => {
            app.show_col_summary = false;
            Task::none()
        }

        Message::AddConsensus => {
            if let Some(arc) = &mut app.document {
                let cons_bytes = helixview_analysis::consensus(
                    arc,
                    helixview_analysis::ConsensusMethod::Plurality,
                );
                let mut aln = (**arc).clone();
                let seq = CoreSequence::new("Consensus", cons_bytes);
                aln.sequences.push(seq);
                *arc = Arc::new(aln);
                app.status = "Consensus row appended".to_string();
            }
            Task::none()
        }

        Message::ToggleFeatures => {
            app.view.show_features = !app.view.show_features;
            Task::none()
        }

        Message::SearchPattern(pattern) => {
            app.view.pattern_matches.clear();
            app.view.pattern_match_idx = 0;
            if let Some(aln) = &app.document {
                if !pattern.is_empty() {
                    let pat_up: Vec<u8> = pattern.bytes()
                        .map(|b| b.to_ascii_uppercase()).collect();
                    for (ri, seq) in aln.sequences.iter().enumerate() {
                        'col: for ci in 0..seq.residues.len().saturating_sub(pat_up.len() - 1) {
                            for (k, &pb) in pat_up.iter().enumerate() {
                                let rb = seq.residues.get(ci + k)
                                    .copied()
                                    .unwrap_or(0)
                                    .to_ascii_uppercase();
                                if !iupac_matches(pb, rb) { continue 'col; }
                            }
                            app.view.pattern_matches.push((ri, ci));
                        }
                    }
                    let n = app.view.pattern_matches.len();
                    if n > 0 {
                        let (row, col) = app.view.pattern_matches[0];
                        app.view.scroll_row = row.min(aln.seq_count().saturating_sub(1));
                        app.view.scroll_col = col;
                        app.status = format!("Pattern '{pattern}': {n} hit(s). Use < > to navigate.");
                    } else {
                        app.status = format!("Pattern '{pattern}': no hits.");
                    }
                }
            }
            Task::none()
        }

        Message::SearchNames(query) => {
            app.view.pattern_matches.clear();
            app.view.pattern_match_idx = 0;
            if let Some(aln) = &app.document {
                if !query.is_empty() {
                    let lower = query.to_lowercase();
                    for (ri, seq) in aln.sequences.iter().enumerate() {
                        if seq.name.to_lowercase().contains(&lower) {
                            app.view.pattern_matches.push((ri, 0));
                        }
                    }
                    let n = app.view.pattern_matches.len();
                    if n > 0 {
                        app.view.scroll_row = app.view.pattern_matches[0].0
                            .min(aln.seq_count().saturating_sub(1));
                        app.status = format!("Name '{query}': {n} match(es). Use < > to navigate.");
                    } else {
                        app.status = format!("Name '{query}': no matches.");
                    }
                }
            }
            Task::none()
        }

        Message::PatternMatchJump(delta) => {
            let n = app.view.pattern_matches.len();
            if n == 0 { return Task::none(); }
            let idx = app.view.pattern_match_idx;
            app.view.pattern_match_idx =
                ((idx as isize + delta).rem_euclid(n as isize)) as usize;
            if let Some(aln) = &app.document {
                let (row, col) = app.view.pattern_matches[app.view.pattern_match_idx];
                app.view.scroll_row = row.min(aln.seq_count().saturating_sub(1));
                app.view.scroll_col = col;
                app.status = format!(
                    "Hit {} / {}",
                    app.view.pattern_match_idx + 1,
                    n,
                );
            }
            Task::none()
        }

        Message::TitleRightClick { idx } => {
            app.title_menu = Some(idx);
            if let Some(aln) = &app.document {
                app.rename_input = aln.sequences
                    .get(idx)
                    .map(|s| s.name.clone())
                    .unwrap_or_default();
            }
            Task::none()
        }

        Message::TitleMenuClose => {
            app.title_menu = None;
            app.rename_input.clear();
            Task::none()
        }

        Message::RenameInput(s) => {
            app.rename_input = s;
            Task::none()
        }

        Message::RenameSequence { idx, name } => {
            if let Some(arc) = &mut app.document {
                let mut aln = (**arc).clone();
                if let Some(seq) = aln.sequences.get_mut(idx) {
                    seq.name = name;
                }
                *arc = Arc::new(aln);
            }
            app.title_menu = None;
            app.rename_input.clear();
            Task::none()
        }

        Message::ToggleLockSequence { idx } => {
            if let Some(arc) = &mut app.document {
                let mut aln = (**arc).clone();
                if let Some(seq) = aln.sequences.get_mut(idx) {
                    seq.locked = !seq.locked;
                    app.status = format!(
                        "Sequence '{}' {}",
                        seq.name,
                        if seq.locked { "locked" } else { "unlocked" }
                    );
                }
                *arc = Arc::new(aln);
            }
            app.title_menu = None;
            Task::none()
        }

        Message::DuplicateSequence { idx } => {
            if let Some(arc) = &mut app.document {
                let mut aln = (**arc).clone();
                if idx < aln.sequences.len() {
                    let mut dup = aln.sequences[idx].clone();
                    dup.name = format!("{} (copy)", dup.name);
                    aln.sequences.insert(idx + 1, dup);
                    *arc = Arc::new(aln);
                    app.status = format!("Duplicated sequence at row {}", idx + 1);
                }
            }
            app.title_menu = None;
            Task::none()
        }

        Message::DeleteSequence { idx } => {
            if let Some(arc) = &mut app.document {
                let mut aln = (**arc).clone();
                if idx < aln.sequences.len() {
                    let name = aln.sequences[idx].name.clone();
                    aln.sequences.remove(idx);
                    *arc = Arc::new(aln);
                    app.view.selected.remove(&idx);
                    app.history = EditHistory::new(100);
                    app.status = format!("Deleted sequence '{name}'");
                }
            }
            app.title_menu = None;
            Task::none()
        }

        Message::ToggleRestrMap => {
            app.show_restr_map = !app.show_restr_map;
            Task::none()
        }

        Message::OpenColorEditor => {
            app.show_color_editor = true;
            Task::none()
        }

        Message::CloseColorEditor => {
            app.show_color_editor = false;
            Task::none()
        }

        Message::SelectColorResidue(b, is_nuc) => {
            app.color_editor_selected = Some((b, is_nuc));
            Task::none()
        }

        Message::SetResidueColor { residue, is_nuc, channel, value } => {
            let mut table = (*app.color_table).clone();
            let cur = if is_nuc { table.nuc_color(residue) } else { table.aa_color(residue) };
            let new_color = match channel {
                0 => iced::Color { r: value, ..cur },
                1 => iced::Color { g: value, ..cur },
                _ => iced::Color { b: value, ..cur },
            };
            if is_nuc { table.set_nuc(residue, new_color); }
            else       { table.set_aa(residue, new_color); }
            app.color_table = Arc::new(table);
            Task::none()
        }

        Message::ResetColorTable => {
            app.color_table = Arc::new(ColorTable::default());
            Task::none()
        }

        Message::ShowConservation => {
            app.show_conservation = true;
            Task::none()
        }

        Message::CloseConservation => {
            app.show_conservation = false;
            Task::none()
        }

        Message::ConservationThreshold(t) => {
            app.conservation_threshold = t.clamp(0.0, 1.0);
            Task::none()
        }

        Message::ConservationMinWidth(w) => {
            app.conservation_min_width = w.max(1);
            Task::none()
        }

        Message::JumpToColumn(col) => {
            app.view.scroll_col = col;
            app.show_mutual_info = false;
            Task::none()
        }

        Message::RunExternalAlign(aligner) => {
            if let Some(aln) = &app.document {
                let aln_clone = (**aln).clone();
                let selected: Vec<usize> = app.view.selected.iter().copied().collect();
                app.status = format!("Running {}…", aligner.display_name());
                return Task::perform(
                    async move {
                        tokio::task::spawn_blocking(move || -> Result<SeqAlignment, String> {
                            helixview_external::run_aligner(&aln_clone, &selected, aligner)
                        })
                        .await
                        .unwrap_or_else(|e| Err(e.to_string()))
                    },
                    Message::ExternalAlignDone,
                );
            }
            Task::none()
        }

        Message::ExternalAlignDone(Ok(new_aln)) => {
            let n = new_aln.seq_count();
            let c = new_aln.col_count();
            app.status = format!("Alignment complete: {n} sequences × {c} columns");
            app.document = Some(Arc::new(new_aln));
            app.view     = ViewState::default();
            app.history  = EditHistory::new(100);
            Task::none()
        }

        Message::ExternalAlignDone(Err(e)) => {
            app.status = format!("Alignment failed: {e}");
            Task::none()
        }

        Message::EditSequenceRaw { idx } => {
            if let Some(aln) = &app.document {
                if let Some(seq) = aln.sequences.get(idx) {
                    let raw = seq.raw_sequence();
                    app.seq_editor_content = Some(text_editor::Content::with_text(&raw));
                    app.seq_editor_idx = Some(idx);
                }
            }
            Task::none()
        }

        Message::SeqEditorAction(action) => {
            if let Some(content) = &mut app.seq_editor_content {
                content.perform(action);
            }
            Task::none()
        }

        Message::SeqEditorCommit => {
            if let (Some(idx), Some(content), Some(arc)) =
                (app.seq_editor_idx, &app.seq_editor_content, &mut app.document)
            {
                let new_raw: String = content.text();
                // Strip whitespace/newlines the user may have introduced.
                let new_bytes: Vec<u8> = new_raw.bytes()
                    .filter(|&b| !b.is_ascii_whitespace())
                    .map(|b| b.to_ascii_uppercase())
                    .collect();
                let mut aln = (**arc).clone();
                // Apply new residues (drop the mutable borrow before borrowing again).
                let seq_name = if let Some(seq) = aln.sequences.get_mut(idx) {
                    seq.residues = new_bytes.clone();
                    seq.name.clone()
                } else {
                    String::new()
                };
                // Pad other sequences to the new maximum column count.
                let new_col_count = aln.sequences.iter().map(|s| s.residues.len()).max().unwrap_or(0);
                for (i, s) in aln.sequences.iter_mut().enumerate() {
                    if i != idx {
                        let pad = new_col_count.saturating_sub(s.residues.len());
                        s.residues.extend(std::iter::repeat(b'-').take(pad));
                    }
                }
                if !seq_name.is_empty() {
                    app.status = format!("Updated '{}': {} residues", seq_name, new_bytes.len());
                }
                *arc = Arc::new(aln);
                app.history = EditHistory::new(100);
            }
            app.seq_editor_idx     = None;
            app.seq_editor_content = None;
            Task::none()
        }

        Message::SeqEditorClose => {
            app.seq_editor_idx     = None;
            app.seq_editor_content = None;
            Task::none()
        }

        Message::ShowDotPlot => {
            // Use the first two selected sequences; fall back to 0 and 1.
            let mut sel: Vec<usize> = app.view.selected.iter().copied().collect();
            sel.sort_unstable();
            if let Some(&a) = sel.first() { app.dot_plot_seq_a = a; }
            if let Some(&b) = sel.get(1)  { app.dot_plot_seq_b = b; }
            // Ensure b != a and indices are in range.
            if let Some(aln) = &app.document {
                let n = aln.seq_count();
                app.dot_plot_seq_a = app.dot_plot_seq_a.min(n.saturating_sub(1));
                app.dot_plot_seq_b = app.dot_plot_seq_b.min(n.saturating_sub(1));
                if app.dot_plot_seq_b == app.dot_plot_seq_a && n > 1 {
                    app.dot_plot_seq_b = if app.dot_plot_seq_a == 0 { 1 } else { 0 };
                }
            }
            app.show_dot_plot = true;
            Task::none()
        }

        Message::CloseDotPlot => {
            app.show_dot_plot = false;
            Task::none()
        }

        Message::DotPlotSeqA(idx) => {
            app.dot_plot_seq_a = idx;
            Task::none()
        }

        Message::DotPlotSeqB(idx) => {
            app.dot_plot_seq_b = idx;
            Task::none()
        }

        Message::DotPlotWindow(w) => {
            app.dot_plot_window    = w.max(1);
            app.dot_plot_min_match = app.dot_plot_min_match.min(app.dot_plot_window);
            Task::none()
        }

        Message::DotPlotMinMatch(m) => {
            app.dot_plot_min_match = m.clamp(1, app.dot_plot_window);
            Task::none()
        }

        Message::SequenceMoved { from, to } => {
            if let Some(arc) = &mut app.document {
                let mut aln = (**arc).clone();
                let cmd = Box::new(MoveSequence { from, to });
                app.history.execute(cmd, &mut aln);
                *arc = Arc::new(aln);
                // Keep the selection on the moved row.
                app.view.selected.clear();
                app.view.selected.insert(to);
                app.view.last_clicked = Some(to);
                app.status = format!("Moved sequence to row {}  (Ctrl+Z to undo)", to + 1);
            }
            Task::none()
        }

        Message::Undo => {
            if let Some(arc) = &mut app.document {
                if app.history.can_undo() {
                    let mut aln = (**arc).clone();
                    app.history.undo(&mut aln);
                    // After undo the command is on the redo stack; its description
                    // is the most recently pushed entry there.
                    let desc = app.history
                        .redo_description()
                        .unwrap_or("edit")
                        .to_owned();
                    *arc = Arc::new(aln);
                    app.status = format!("Undo: {desc}");
                }
            }
            Task::none()
        }

        Message::Redo => {
            if let Some(arc) = &mut app.document {
                if app.history.can_redo() {
                    let mut aln = (**arc).clone();
                    app.history.redo(&mut aln);
                    let desc = app.history
                        .undo_description()
                        .unwrap_or("edit")
                        .to_owned();
                    *arc = Arc::new(aln);
                    app.status = format!("Redo: {desc}");
                }
            }
            Task::none()
        }

        // ── NCBI Entrez fetch ────────────────────────────────────────────────
        Message::ToggleFetchBar => {
            app.view.show_fetch_bar = !app.view.show_fetch_bar;
            if !app.view.show_fetch_bar {
                app.view.fetch_input.clear();
            }
            Task::none()
        }

        Message::FetchInput(s) => {
            app.view.fetch_input = s;
            Task::none()
        }

        Message::FetchAccession => {
            let acc = app.view.fetch_input.trim().to_string();
            if acc.is_empty() { return Task::none(); }
            app.status = format!("Fetching '{acc}' from NCBI…");
            Task::perform(
                entrez::fetch_accession(acc),
                Message::EntrezResult,
            )
        }

        Message::EntrezResult(result) => {
            match result {
                Ok(fetched) => {
                    let n = fetched.seq_count();
                    let names: Vec<_> = fetched.sequences.iter()
                        .map(|s| s.name.as_str())
                        .take(3)
                        .collect();
                    let label = names.join(", ");
                    match &mut app.document {
                        None => {
                            app.document = Some(Arc::new(fetched));
                            app.history  = EditHistory::new(100);
                            app.view     = ViewState::default();
                        }
                        Some(arc) => {
                            let mut aln = (**arc).clone();
                            for seq in fetched.sequences {
                                aln.push(seq);
                            }
                            *arc = Arc::new(aln);
                        }
                    }
                    app.view.show_fetch_bar = false;
                    app.view.fetch_input.clear();
                    app.status = format!("Fetched {n} sequence(s): {label}");
                }
                Err(e) => {
                    app.status = format!("Fetch error: {e}");
                }
            }
            Task::none()
        }

        // ── Export ───────────────────────────────────────────────────────────
        Message::ExportFastaDialog => {
            Task::perform(
                async { rfd::AsyncFileDialog::new()
                    .add_filter("FASTA", &["fa", "fasta", "fna", "faa"])
                    .save_file().await
                    .map(|h| h.path().to_path_buf())
                },
                |p| Message::ExportFastaChosen { path: p, raw: false },
            )
        }

        Message::ExportRawFastaDialog => {
            Task::perform(
                async { rfd::AsyncFileDialog::new()
                    .add_filter("FASTA", &["fa", "fasta", "fna", "faa"])
                    .save_file().await
                    .map(|h| h.path().to_path_buf())
                },
                |p| Message::ExportFastaChosen { path: p, raw: true },
            )
        }

        Message::ExportFastaChosen { path: None, .. } => Task::none(),

        Message::ExportFastaChosen { path: Some(path), raw } => {
            if let Some(arc) = &app.document {
                // Collect sequences to export: selected rows, or all.
                let indices: Vec<usize> = if app.view.selected.is_empty() {
                    (0..arc.seq_count()).collect()
                } else {
                    app.view.selected.iter().copied().collect()
                };

                let mut export_aln = SeqAlignment::new(&arc.name);
                for idx in &indices {
                    if let Some(seq) = arc.sequences.get(*idx) {
                        let mut s = seq.clone();
                        if raw {
                            s.residues.retain(|&b| !helixview_core::sequence::is_gap(b));
                            s.gap_locks.clear();
                        }
                        export_aln.push(s);
                    }
                }

                match helixview_formats::fasta::write_file(&export_aln, &path) {
                    Ok(()) => app.status = format!(
                        "Exported {} sequence(s) to {}{}",
                        indices.len(),
                        path.display(),
                        if raw { " (gaps stripped)" } else { "" }
                    ),
                    Err(e) => app.status = format!("Export error: {e}"),
                }
            }
            Task::none()
        }

        Message::ExportIdentityCsv => {
            Task::perform(
                async { rfd::AsyncFileDialog::new()
                    .add_filter("CSV", &["csv"])
                    .save_file().await
                    .map(|h| h.path().to_path_buf())
                },
                Message::ExportCsvChosen,
            )
        }

        Message::ExportCsvChosen(None) => Task::none(),

        Message::ExportCsvChosen(Some(path)) => {
            if let Some(arc) = &app.document {
                let n = arc.seq_count().min(200);
                let mut csv = String::new();
                // Header row
                csv.push(',');
                for j in 0..n {
                    csv.push_str(&format!("{}", arc.sequences[j].name.replace(',', "_")));
                    if j + 1 < n { csv.push(','); }
                }
                csv.push('\n');
                for i in 0..n {
                    csv.push_str(&arc.sequences[i].name.replace(',', "_"));
                    for j in 0..n {
                        let pct = if i == j {
                            100.0f64
                        } else {
                            let a: Vec<u8> = arc.sequences[i].residues.iter()
                                .filter(|&&b| !helixview_core::sequence::is_gap(b))
                                .copied().collect();
                            let b: Vec<u8> = arc.sequences[j].residues.iter()
                                .filter(|&&b| !helixview_core::sequence::is_gap(b))
                                .copied().collect();
                            helixview_analysis::pairwise_identity(&a, &b) * 100.0
                        };
                        csv.push_str(&format!(",{:.1}", pct));
                    }
                    csv.push('\n');
                }
                match std::fs::write(&path, &csv) {
                    Ok(()) => app.status = format!("Identity matrix exported to {}", path.display()),
                    Err(e) => app.status = format!("CSV export error: {e}"),
                }
            }
            Task::none()
        }

        // ── Feature annotation ───────────────────────────────────────────────
        Message::ToggleAnnotateBar => {
            app.view.show_annotate_bar = !app.view.show_annotate_bar;
            Task::none()
        }

        Message::AnnotateName(s) => {
            app.view.annotate_name = s;
            Task::none()
        }

        Message::AnnotateType(s) => {
            app.view.annotate_type = s;
            Task::none()
        }

        Message::ApplyAnnotation => {
            if app.view.selected_cols.is_empty() || app.view.annotate_name.trim().is_empty() {
                app.status = "Select columns and enter a feature name first.".to_string();
                return Task::none();
            }
            if let Some(arc) = &mut app.document {
                let col_lo = *app.view.selected_cols.iter().next().unwrap();
                let col_hi = *app.view.selected_cols.iter().next_back().unwrap();
                let feat_type = helixview_core::Feature::type_from_str(&app.view.annotate_type);

                let target_rows: Vec<usize> = if app.view.selected.is_empty() {
                    (0..arc.seq_count()).collect()
                } else {
                    app.view.selected.iter().copied().collect()
                };

                let mut aln = (**arc).clone();
                let mut annotated = 0usize;
                for idx in &target_rows {
                    if let Some(seq) = aln.sequences.get_mut(*idx) {
                        // Count non-gap residues before each column boundary.
                        let true_start = seq.residues[..col_lo.min(seq.residues.len())]
                            .iter().filter(|&&b| !helixview_core::sequence::is_gap(b)).count();
                        let true_end = seq.residues[..=col_hi.min(seq.residues.len().saturating_sub(1))]
                            .iter().filter(|&&b| !helixview_core::sequence::is_gap(b)).count()
                            .saturating_sub(1);

                        let mut feat = helixview_core::Feature::new(
                            app.view.annotate_name.trim(),
                            true_start,
                            true_end,
                        );
                        feat.feature_type = feat_type.clone();
                        feat.color = feat_type.default_color();
                        seq.features.push(feat);
                        annotated += 1;
                    }
                }
                *arc = Arc::new(aln);
                app.view.show_annotate_bar = false;
                app.view.annotate_name.clear();
                app.status = format!(
                    "Annotated {} sequence(s): '{}' cols {}–{}",
                    annotated,
                    app.view.annotate_type,
                    col_lo + 1, col_hi + 1
                );
            }
            Task::none()
        }

        // ── Protein view / translation toggle ────────────────────────────────
        Message::SetProteinViewFrame(frame) => {
            if let Some(arc) = &app.document {
                let protein_aln = make_protein_alignment(arc, frame);
                app.protein_view_doc   = Some(Arc::new(protein_aln));
                app.protein_view_frame = frame;
                app.view.protein_view_frame = Some(frame);
                let label: i8 = if frame < 3 { (frame as i8) + 1 } else { -((frame as i8) - 2) };
                app.status = format!("Protein view — frame {:+}", label);
            }
            Task::none()
        }

        Message::ExitProteinView => {
            app.protein_view_doc = None;
            app.view.protein_view_frame = None;
            app.status = "Returned to DNA view.".to_string();
            Task::none()
        }

        // ── BLAST ─────────────────────────────────────────────────────────────
        Message::ToggleBlastPanel => {
            app.view.show_blast_panel = !app.view.show_blast_panel;
            Task::none()
        }

        Message::BlastProgramChanged(p) => { app.view.blast_program = p; Task::none() }
        Message::BlastDatabaseChanged(d) => { app.view.blast_database = d; Task::none() }

        Message::SubmitBlast => {
            if let Some(arc) = &app.document {
                let idx = app.view.selected.iter().next().copied()
                    .or(app.view.hover_row)
                    .unwrap_or(0);
                if let Some(seq) = arc.sequences.get(idx) {
                    let residues = seq.residues.iter()
                        .filter(|&&b| !helixview_core::sequence::is_gap(b))
                        .copied().collect::<Vec<u8>>();
                    let program  = app.view.blast_program.clone();
                    let database = app.view.blast_database.clone();
                    app.status = format!("Submitting {} to NCBI BLAST…", seq.name);
                    app.blast  = BlastState::Idle;
                    return Task::perform(
                        async move { blast::submit(&residues, &program, &database).await },
                        Message::BlastSubmitted,
                    );
                }
            }
            Task::none()
        }

        Message::BlastSubmitted(Ok(rid)) => {
            app.status = format!("BLAST submitted — RID: {rid}  (click 'Check Results' to poll)");
            app.blast  = BlastState::Submitted(rid);
            app.view.show_blast_panel = false;
            Task::none()
        }

        Message::BlastSubmitted(Err(e)) => {
            app.status = format!("BLAST submit failed: {e}");
            app.blast  = BlastState::Failed(e);
            Task::none()
        }

        Message::BlastCheck => {
            if let BlastState::Submitted(ref rid) = app.blast {
                let rid = rid.clone();
                app.status = "Polling NCBI BLAST…".to_string();
                return Task::perform(
                    async move { blast::poll(&rid).await },
                    Message::BlastPollResult,
                );
            }
            Task::none()
        }

        Message::BlastPollResult(Ok(true)) => {
            if let BlastState::Submitted(ref rid) = app.blast.clone() {
                let rid = rid.clone();
                app.status = "BLAST results ready — fetching…".to_string();
                return Task::perform(
                    async move { blast::fetch_hits(&rid, 25).await },
                    Message::BlastResultsFetched,
                );
            }
            Task::none()
        }

        Message::BlastPollResult(Ok(false)) => {
            app.status = "BLAST still running — try again in a few seconds.".to_string();
            Task::none()
        }

        Message::BlastPollResult(Err(e)) => {
            app.blast  = BlastState::Failed(e.clone());
            app.status = format!("BLAST poll error: {e}");
            Task::none()
        }

        Message::BlastResultsFetched(Ok(hits)) => {
            let n = hits.len();
            app.status = format!("BLAST complete — {n} hit(s) found.");
            app.blast  = BlastState::Complete(hits);
            Task::none()
        }

        Message::BlastResultsFetched(Err(e)) => {
            app.blast  = BlastState::Failed(e.clone());
            app.status = format!("BLAST fetch error: {e}");
            Task::none()
        }

        Message::BlastImportHit(accession) => {
            app.status = format!("Fetching hit '{accession}' from NCBI…");
            Task::perform(
                entrez::fetch_accession(accession),
                Message::EntrezResult,
            )
        }

        Message::CloseBlast => {
            app.blast = BlastState::Idle;
            Task::none()
        }

        // ── Toolbar ribbon ────────────────────────────────────────────────────
        Message::SetToolbarTab(tab) => {
            app.view.toolbar_tab = tab;
            Task::none()
        }

        // ── ABI trace viewer ──────────────────────────────────────────────────
        Message::OpenTraceDialog => {
            Task::perform(pick_trace_file(), Message::TraceFileChosen)
        }

        Message::TraceFileChosen(None) => Task::none(),

        Message::TraceFileChosen(Some(path)) => {
            app.status = format!("Loading ABI trace: {}…", path.display());
            Task::perform(
                async move {
                    let data = tokio::fs::read(&path).await
                        .map_err(|e| e.to_string())?;
                    helixview_formats::abi::parse(&data)
                },
                Message::TraceLoaded,
            )
        }

        Message::TraceLoaded(Ok(trace)) => {
            app.status = format!(
                "Loaded ABI trace: {} bases, {} samples",
                trace.bases.len(), trace.num_samples,
            );
            app.trace_scroll = 0.0;
            app.trace_zoom   = 1.0;
            // Choose default quant channels from FWO_ order.
            app.trace_quant_ch_a = trace.channel_bases[0];
            app.trace_quant_ch_b = trace.channel_bases[1];
            app.trace_data = Some(Arc::new(trace));
            Task::none()
        }

        Message::TraceLoaded(Err(e)) => {
            app.status = format!("Failed to load ABI trace: {e}");
            Task::none()
        }

        Message::CloseTrace => {
            app.trace_data = None;
            Task::none()
        }

        Message::TraceScroll(v) => {
            app.trace_scroll = v.clamp(0.0, 1.0);
            Task::none()
        }

        Message::TraceZoom(v) => {
            app.trace_zoom = v.clamp(1.0, 16.0);
            Task::none()
        }

        Message::TraceQuantChA(b) => {
            app.trace_quant_ch_a = b;
            Task::none()
        }

        Message::TraceQuantChB(b) => {
            app.trace_quant_ch_b = b;
            Task::none()
        }

        // ── Tree viewer ───────────────────────────────────────────────────────

        Message::BuildTreeFromAlignment => {
            if let Some(aln) = &app.document {
                let n = aln.seq_count();
                if n < 2 {
                    app.status = "Need at least 2 sequences to build a tree.".to_string();
                    return Task::none();
                }
                app.computing = Some(format!("Building NJ tree ({n} seqs)…"));
                app.status = format!("Building neighbour-joining tree for {n} sequences…");
                let aln = Arc::clone(aln);
                Task::perform(
                    async move {
                        tokio::task::spawn_blocking(move || {
                            helixview_analysis::nj_tree_from_alignment(&aln)
                        })
                        .await
                        .map(|t| vec![t])
                        .map_err(|e| e.to_string())
                    },
                    Message::TreeLoaded,
                )
            } else {
                app.status = "Open an alignment first to build a tree.".to_string();
                Task::none()
            }
        }

        Message::OpenTreeDialog => {
            Task::perform(pick_tree_file(), Message::TreeFileChosen)
        }

        Message::TreeFileChosen(None) => Task::none(),

        Message::TreeFileChosen(Some(path)) => {
            app.status = format!("Loading tree {}…", path.display());
            Task::perform(
                async move {
                    helixview_formats::newick::parse_file(&path)
                        .map_err(|e| e.to_string())
                },
                Message::TreeLoaded,
            )
        }

        Message::TreeLoaded(Err(e)) => {
            app.status = format!("Tree error: {e}");
            Task::none()
        }

        Message::TreeLoaded(Ok(trees)) => {
            app.computing = None;
            if let Some(tree) = trees.into_iter().next() {
                let lc = tree.leaf_count();
                app.status = format!("Tree ready: {} leaves", lc);
                app.tree_data = Some(Arc::new(tree));
                app.tree_selected.clear();
            }
            Task::none()
        }

        Message::CloseTree => {
            app.tree_data = None;
            app.tree_selected.clear();
            Task::none()
        }

        Message::TreeNodeClicked(name) => {
            if app.tree_selected.contains(&name) {
                app.tree_selected.remove(&name);
            } else {
                app.tree_selected.insert(name);
            }
            Task::none()
        }

        Message::TreeClearSelection => {
            app.tree_selected.clear();
            Task::none()
        }

        // ── Composition & Tm ─────────────────────────────────────────────────

        Message::ShowComposition  => { app.show_composition = true;  Task::none() }
        Message::CloseComposition => { app.show_composition = false; Task::none() }
        Message::ShowOligoTm      => { app.show_oligo_tm = true;  Task::none() }
        Message::CloseOligoTm     => { app.show_oligo_tm = false; Task::none() }

        // ── Hydrophobicity profile ────────────────────────────────────────────
        Message::ShowHydrophobicity  => { app.show_hydrophobicity = true;  Task::none() }
        Message::CloseHydrophobicity => { app.show_hydrophobicity = false; Task::none() }
        Message::HydroWindow(w)      => { app.hydro_window = w.max(3).min(25); Task::none() }

        // ── PDF export ────────────────────────────────────────────────────────
        Message::ExportPdfDialog => {
            Task::perform(save_pdf_dialog(), Message::PdfExportPath)
        }
        Message::PdfExportPath(None) => Task::none(),
        Message::PdfExportPath(Some(path)) => {
            if let Some(aln) = &app.document {
                let aln = Arc::clone(aln);
                let opts = app.shaded_export_opts.clone();
                let col_end = aln.col_count();
                let color_table = Arc::clone(&app.color_table);
                Task::perform(
                    async move {
                        export_pdf_task(aln, opts, 0, col_end, color_table, path).await
                    },
                    Message::PdfExportDone,
                )
            } else {
                Task::none()
            }
        }
        Message::PdfExportDone(Ok(())) => {
            app.status = "PDF exported.".to_string();
            Task::none()
        }
        Message::PdfExportDone(Err(e)) => {
            app.status = format!("PDF export failed: {e}");
            Task::none()
        }

        // ── Text export ───────────────────────────────────────────────────────
        Message::ShowTextExport    => { app.show_text_export = true;  Task::none() }
        Message::CloseTextExport   => { app.show_text_export = false; Task::none() }
        Message::TextExportResPerRow(v) => { app.text_export_opts.residues_per_row = v; Task::none() }
        Message::TextExportTitleChars(v) => { app.text_export_opts.title_chars = v; Task::none() }
        Message::TextExportShowRuler(v) => { app.text_export_opts.show_ruler = v; Task::none() }
        Message::TextExportNumberLines(v) => { app.text_export_opts.number_lines = v; Task::none() }
        Message::TextExportSave => {
            Task::perform(save_text_dialog(), Message::TextExportPath)
        }
        Message::TextExportPath(None) => Task::none(),
        Message::TextExportPath(Some(path)) => {
            if let Some(aln) = &app.document {
                let aln = Arc::clone(aln);
                let opts = app.text_export_opts.clone();
                Task::perform(
                    async move {
                        let text = crate::text_export::format_alignment(&aln, &opts);
                        tokio::task::spawn_blocking(move || {
                            std::fs::write(&path, text).map_err(|e| e.to_string())
                        })
                        .await
                        .map_err(|e| e.to_string())?
                    },
                    Message::TextExportDone,
                )
            } else {
                Task::none()
            }
        }
        Message::TextExportDone(Ok(())) => {
            app.status = "Alignment text saved.".to_string();
            Task::none()
        }
        Message::TextExportDone(Err(e)) => {
            app.status = format!("Text export failed: {e}");
            Task::none()
        }

        // ── Mutual information ────────────────────────────────────────────────
        Message::ShowMutualInfo  => {
            app.show_mutual_info = true;
            let stale = app.mi_aln_ptr.map_or(true, |p| {
                app.document.as_ref().map_or(true, |a| Arc::as_ptr(a) as usize != p)
            });
            if app.document.is_some() && (app.mi_result.is_none() || stale) && !app.mi_running {
                Task::done(Message::RunMutualInfo)
            } else {
                Task::none()
            }
        }
        Message::CloseMutualInfo => { app.show_mutual_info = false; Task::none() }
        Message::MiMinObs(v)    => { app.mi_min_obs = v; Task::none() }
        Message::MiTopN(v)      => { app.mi_top_n   = v; Task::none() }
        Message::MiStemMinWcTypes(v) => { app.mi_stem_min_wc_types = v; Task::none() }
        Message::MiStemMinWcFrac(v)  => { app.mi_stem_min_wc_frac  = v; Task::none() }
        Message::MiSetTab(tab)      => { app.mi_tab = tab; Task::none() }
        Message::TogglePairingArcs  => { app.show_pairing_arcs = !app.show_pairing_arcs; Task::none() }
        Message::PairingArcThreshold(v) => { app.pairing_arc_threshold = v; Task::none() }
        Message::RunMutualInfo  => {
            if let Some(aln) = &app.document {
                if aln.seq_count() < 2 {
                    app.status = "MI requires at least 2 sequences — load an alignment first.".to_string();
                    return Task::none();
                }
                app.mi_running = true;
                app.computing = Some("Computing mutual information…".to_string());
                let aln = Arc::clone(aln);
                let min_obs = app.mi_min_obs;
                let top_n   = app.mi_top_n;
                Task::perform(
                    async move {
                        tokio::task::spawn_blocking(move || {
                            helixview_analysis::mutual_information(&aln, min_obs, top_n)
                        })
                        .await
                        .map(Arc::new)
                        .map_err(|e| e.to_string())
                    },
                    Message::MiDone,
                )
            } else {
                Task::none()
            }
        }
        Message::MiDone(Ok(result)) => {
            app.mi_running = false;
            app.computing = None;
            let n_seqs = app.document.as_ref().map_or(0, |a| a.seq_count());
            app.status = format!(
                "MI done — {} seqs · {} cols · {} top pairs · max MI {:.3}",
                n_seqs, result.n_cols, result.top_pairs.len(), result.max_mi()
            );
            app.mi_result  = Some(result);
            app.mi_aln_ptr = app.document.as_ref().map(|a| Arc::as_ptr(a) as usize);
            Task::none()
        }
        Message::MiDone(Err(e)) => {
            app.mi_running = false;
            app.computing = None;
            app.status = format!("MI failed: {e}");
            Task::none()
        }

        // ── Async identity matrix ─────────────────────────────────────────────
        Message::ComputeIdentityMatrix => {
            if let Some(aln) = &app.document {
                app.computing = Some("Computing identity matrix…".to_string());
                app.identity_matrix_data = None;
                let aln = Arc::clone(aln);
                Task::perform(
                    async move {
                        tokio::task::spawn_blocking(move || {
                            let n = aln.seq_count().min(200);
                            let matrix: Vec<Vec<f64>> = (0..n)
                                .map(|i| {
                                    let a: Vec<u8> = aln.sequences[i].residues.iter()
                                        .filter(|&&b| !matches!(b, b'-'|b'~'|b'.'))
                                        .copied().collect();
                                    (0..n).map(|j| {
                                        if i == j { return 1.0; }
                                        let b_seq: Vec<u8> = aln.sequences[j].residues.iter()
                                            .filter(|&&b| !matches!(b, b'-'|b'~'|b'.'))
                                            .copied().collect();
                                        helixview_analysis::pairwise_identity(&a, &b_seq)
                                    }).collect()
                                })
                                .collect();
                            Arc::new(matrix)
                        })
                        .await
                        .unwrap_or_else(|_| Arc::new(vec![]))
                    },
                    Message::IdentityMatrixDone,
                )
            } else {
                Task::none()
            }
        }
        Message::IdentityMatrixDone(data) => {
            app.computing = None;
            app.identity_matrix_data = Some(data);
            Task::none()
        }

        // ── Async pairwise ────────────────────────────────────────────────────
        Message::PairwiseDone(Ok(result)) => {
            app.computing = None;
            app.status = format!(
                "Pairwise: score={}, identity={:.1}%",
                result.score, result.identity * 100.0
            );
            app.pairwise_result = Some(result);
            Task::none()
        }
        Message::PairwiseDone(Err(e)) => {
            app.computing = None;
            app.status = format!("Pairwise failed: {e}");
            Task::none()
        }

        // ── Local BLAST DB ────────────────────────────────────────────────────
        Message::CreateBlastDbDialog => {
            Task::perform(
                async {
                    rfd::AsyncFileDialog::new()
                        .set_title("Select FASTA file to build BLAST database…")
                        .add_filter("FASTA", &["fasta", "fa", "fna", "faa"])
                        .add_filter("All files", &["*"])
                        .pick_file()
                        .await
                        .map(|h| h.path().to_path_buf())
                },
                Message::CreateBlastDbFasta,
            )
        }
        Message::CreateBlastDbFasta(None) => Task::none(),
        Message::CreateBlastDbFasta(Some(fasta_path)) => {
            app.computing = Some("Creating BLAST database…".to_string());
            let db_path = fasta_path.with_extension("");
            // Infer dbtype: check extension first, then sniff content.
            let ext = fasta_path.extension()
                .and_then(|e| e.to_str())
                .unwrap_or("")
                .to_ascii_lowercase();
            let ext_is_protein = matches!(ext.as_str(), "faa" | "pep" | "prot" | "aa");
            // Content-sniff: if any residue in the first sequence is protein-only
            // (E, F, I, L, P, Q, Z — not present in nucleotide alphabets), it's protein.
            let content_is_protein = std::fs::read(&fasta_path).ok().map_or(false, |bytes| {
                let protein_only: &[u8] = b"EFILPQZefilpqz";
                bytes.iter().any(|b| protein_only.contains(b))
            });
            let dbtype = if ext_is_protein || content_is_protein { "prot" } else { "nucl" };
            let dbtype = dbtype.to_string();
            Task::perform(
                async move {
                    tokio::task::spawn_blocking(move || {
                        let out = std::process::Command::new("makeblastdb")
                            .args([
                                "-in",  fasta_path.to_str().unwrap_or(""),
                                "-dbtype", &dbtype,
                                "-out", db_path.to_str().unwrap_or(""),
                                "-parse_seqids",
                            ])
                            .output()
                            .map_err(|e| format!("makeblastdb not found: {e}"))?;
                        if out.status.success() {
                            let stdout = String::from_utf8_lossy(&out.stdout);
                            Ok(format!("BLAST DB created. {}", stdout.lines().next().unwrap_or("")))
                        } else {
                            Err(String::from_utf8_lossy(&out.stderr).to_string())
                        }
                    })
                    .await
                    .map_err(|e| e.to_string())
                    .and_then(|r| r)
                },
                Message::CreateBlastDbDone,
            )
        }
        Message::CreateBlastDbDone(Ok(msg)) => {
            app.computing = None;
            app.status = msg;
            Task::none()
        }
        Message::CreateBlastDbDone(Err(e)) => {
            app.computing = None;
            app.status = format!("BLAST DB failed: {e}");
            Task::none()
        }

        // ── Accessory apps ────────────────────────────────────────────────────
        Message::ShowAccessories   => { app.show_accessories = true;  Task::none() }
        Message::CloseAccessories  => {
            app.show_accessories = false;
            app.editing_accessory = None;
            app.acc_edit_name.clear();
            app.acc_edit_cmd.clear();
            Task::none()
        }
        Message::NewAccessory => {
            app.editing_accessory = None;
            app.acc_edit_name = String::new();
            app.acc_edit_cmd  = String::new();
            Task::none()
        }
        Message::EditAccessory(i) => {
            if let Some(acc) = app.accessories.get(i) {
                app.editing_accessory = Some(i);
                app.acc_edit_name = acc.name.clone();
                app.acc_edit_cmd  = acc.command.clone();
            }
            Task::none()
        }
        Message::DeleteAccessory(i) => {
            if i < app.accessories.len() {
                app.accessories.remove(i);
                let _ = crate::accessories::save_accessories(&app.accessories);
            }
            Task::none()
        }
        Message::SaveAccessory => {
            let acc = AccessoryApp {
                name:    app.acc_edit_name.trim().to_string(),
                command: app.acc_edit_cmd.trim().to_string(),
            };
            if !acc.name.is_empty() && !acc.command.is_empty() {
                if let Some(idx) = app.editing_accessory {
                    if idx < app.accessories.len() {
                        app.accessories[idx] = acc;
                    }
                } else {
                    app.accessories.push(acc);
                }
                let _ = crate::accessories::save_accessories(&app.accessories);
                app.editing_accessory = None;
                app.acc_edit_name.clear();
                app.acc_edit_cmd.clear();
            }
            Task::none()
        }
        Message::AccessoryFieldName(s) => { app.acc_edit_name = s; Task::none() }
        Message::AccessoryFieldCmd(s)  => { app.acc_edit_cmd  = s; Task::none() }
        Message::RunAccessory(i) => {
            let acc = match app.accessories.get(i) {
                Some(a) => a.clone(),
                None    => return Task::none(),
            };
            let seqs_to_export = if let Some(aln) = &app.document {
                // Export selected or all sequences as FASTA.
                let idxs: Vec<usize> = if app.view.selected.is_empty() {
                    (0..aln.seq_count()).collect()
                } else {
                    app.view.selected.iter().copied().collect()
                };
                let mut fasta = String::new();
                for idx in idxs {
                    if let Some(s) = aln.sequences.get(idx) {
                        fasta.push('>');
                        fasta.push_str(&s.name);
                        fasta.push('\n');
                        fasta.push_str(&String::from_utf8_lossy(&s.residues));
                        fasta.push('\n');
                    }
                }
                fasta
            } else {
                String::new()
            };
            app.computing = Some(format!("Running {}…", acc.name));
            Task::perform(
                async move {
                    tokio::task::spawn_blocking(move || {
                        // Write FASTA to temp file.
                        let tmp = std::env::temp_dir().join("helixview_acc_input.fasta");
                        std::fs::write(&tmp, &seqs_to_export)
                            .map_err(|e| e.to_string())?;
                        let cmd_str = acc.command
                            .replace("{input}", tmp.to_str().unwrap_or(""))
                            .replace("{output}", std::env::temp_dir()
                                .join("helixview_acc_output.txt")
                                .to_str().unwrap_or(""));
                        // Run via shell.
                        let out = std::process::Command::new("sh")
                            .args(["-c", &cmd_str])
                            .output()
                            .map_err(|e| e.to_string())?;
                        let stdout = String::from_utf8_lossy(&out.stdout).to_string();
                        let stderr = String::from_utf8_lossy(&out.stderr).to_string();
                        if out.status.success() {
                            Ok(format!("{} done. {}", acc.name,
                                stdout.lines().next().unwrap_or("").to_string()))
                        } else {
                            Err(stderr.lines().next().unwrap_or("Unknown error").to_string())
                        }
                    })
                    .await
                    .map_err(|e| e.to_string())
                    .and_then(|r| r)
                },
                Message::AccessoryDone,
            )
        }
        Message::AccessoryDone(Ok(msg)) => {
            app.computing = None;
            app.status = msg;
            Task::none()
        }
        Message::AccessoryDone(Err(e)) => {
            app.computing = None;
            app.status = format!("Accessory failed: {e}");
            Task::none()
        }

        // ── Command palette ───────────────────────────────────────────────────
        Message::OpenCommandPalette => {
            app.show_command_palette = true;
            app.palette_query.clear();
            app.palette_selected = 0;
            text_input::focus(text_input::Id::new("palette_input"))
        }
        Message::CloseCommandPalette => {
            app.show_command_palette = false;
            Task::none()
        }
        Message::PaletteQuery(q) => {
            app.palette_selected = 0;
            app.palette_query = q;
            Task::none()
        }
        Message::PaletteUp => {
            if app.palette_selected > 0 { app.palette_selected -= 1; }
            Task::none()
        }
        Message::PaletteDown => {
            let n = crate::views::palette_count(&app.palette_query);
            if n > 0 && app.palette_selected + 1 < n {
                app.palette_selected += 1;
            }
            Task::none()
        }
        Message::PaletteSelect(orig_idx) => {
            app.show_command_palette = false;
            if let Some(msg) = crate::views::palette_fire(orig_idx) {
                return update(app, msg);
            }
            Task::none()
        }
        Message::PaletteConfirm => {
            if let Some(orig) = crate::views::palette_resolve(&app.palette_query, app.palette_selected) {
                app.show_command_palette = false;
                if let Some(msg) = crate::views::palette_fire(orig) {
                    return update(app, msg);
                }
            }
            Task::none()
        }

        // ── Help ──────────────────────────────────────────────────────────────
        Message::ShowHelp  => { app.show_help = true;  Task::none() }
        Message::CloseHelp => { app.show_help = false; Task::none() }

        // ── Project file ──────────────────────────────────────────────────────
        Message::SaveProjectDialog => {
            Task::perform(save_project_dialog(), Message::SaveProjectChosen)
        }
        Message::SaveProjectChosen(None) => Task::none(),
        Message::SaveProjectChosen(Some(path)) => {
            let result = crate::project::save(app, &path)
                .map_err(|e| e);
            Task::done(Message::SaveProjectDone(result))
        }
        Message::SaveProjectDone(Ok(())) => {
            app.status = "Project saved.".to_string();
            auto_save_session(app);
            Task::none()
        }
        Message::SaveProjectDone(Err(e)) => {
            app.status = format!("Save project failed: {e}");
            Task::none()
        }
        Message::OpenProjectDialog => {
            Task::perform(open_project_dialog(), Message::OpenProjectChosen)
        }
        Message::OpenProjectChosen(None) => Task::none(),
        Message::OpenProjectChosen(Some(path)) => {
            let result = crate::project::load(app, &path);
            Task::done(Message::OpenProjectDone(result))
        }
        Message::OpenProjectDone(Ok(())) => {
            app.status = "Project loaded.".to_string();
            auto_save_session(app);
            Task::none()
        }
        Message::OpenProjectDone(Err(e)) => {
            app.status = format!("Load project failed: {e}");
            Task::none()
        }

        // ── Preferences ───────────────────────────────────────────────────────
        Message::OpenPrefs   => { app.show_prefs = true;  Task::none() }
        Message::ClosePrefs  => { app.show_prefs = false; Task::none() }
        Message::SavePrefs   => {
            if let Err(e) = crate::prefs::save(&app.prefs) {
                app.status = format!("Could not save preferences: {e}");
            } else {
                app.status = "Preferences saved.".to_string();
            }
            app.show_prefs = false;
            Task::none()
        }
        Message::PrefsUndoDepth(v)         => { app.prefs.undo_depth         = v; Task::none() }
        Message::PrefsDefaultZoom(v)       => { app.prefs.default_zoom       = v; Task::none() }
        Message::PrefsFontSize(v)          => { app.prefs.font_size          = v; Task::none() }
        Message::PrefsShowFeatures(v)      => { app.prefs.show_features      = v; Task::none() }
        Message::PrefsShowRestrMap(v)      => { app.prefs.show_restr_map     = v; Task::none() }
        Message::PrefsIdentityThreshold(v) => { app.prefs.identity_threshold = v; Task::none() }
        Message::PrefsRecentFilesMax(v)    => { app.prefs.recent_files_max   = v; Task::none() }
        Message::PrefsClearRecent          => { app.prefs.recent_files.clear(); Task::none() }
        Message::OpenRecentFile(path)      => {
            app.show_prefs = false;
            app.status = format!("Loading {}…", path.display());
            Task::perform(load_file(path), Message::FileLoaded)
        }

        // ── Shaded graphic export ─────────────────────────────────────────────
        Message::ShowShadedExportDialog  => { app.show_shaded_export = true;  Task::none() }
        Message::CloseShadedExportDialog => { app.show_shaded_export = false; Task::none() }
        Message::ShadedExportSetScheme(s) => {
            app.shaded_export_opts.scheme = s;
            Task::none()
        }
        Message::ShadedExportThreshold(t) => {
            app.shaded_export_opts.threshold = t;
            Task::none()
        }
        Message::ShadedExportToggleRuler(v) => {
            app.shaded_export_opts.show_ruler = v;
            Task::none()
        }
        Message::ShadedExportToggleConsensus(v) => {
            app.shaded_export_opts.show_consensus = v;
            Task::none()
        }
        Message::ShadedExportRun => {
            if app.document.is_some() {
                Task::perform(save_shaded_svg_dialog(), Message::ShadedExportChosen)
            } else {
                Task::none()
            }
        }
        Message::ShadedExportChosen(None) => {
            app.show_shaded_export = false;
            Task::none()
        }
        Message::ShadedExportChosen(Some(path)) => {
            if let Some(aln) = &app.document {
                let aln         = Arc::clone(aln);
                let rows        = app.view.selected.clone();
                let col_start   = app.view.selected_cols.iter().next().copied()
                    .unwrap_or(0);
                let col_end     = app.view.selected_cols.iter().last().copied()
                    .map(|c| c + 1)
                    .unwrap_or_else(|| aln.col_count());
                let color_table = Arc::clone(&app.color_table);
                let opts        = app.shaded_export_opts.clone();
                app.show_shaded_export = false;
                app.status = "Exporting shaded graphic…".to_string();
                Task::perform(
                    async move {
                        let svg = crate::export_shaded::shaded_alignment_svg(
                            &aln, &rows, col_start, col_end,
                            &color_table, &opts.scheme, opts.threshold,
                            opts.show_ruler, opts.show_consensus,
                        );
                        std::fs::write(&path, svg.as_bytes()).map_err(|e| e.to_string())
                    },
                    Message::ShadedExportDone,
                )
            } else {
                Task::none()
            }
        }
        Message::ShadedExportDone(Ok(())) => {
            app.status = "Shaded graphic exported.".to_string();
            Task::none()
        }
        Message::ShadedExportDone(Err(e)) => {
            app.status = format!("Shaded export failed: {e}");
            Task::none()
        }

        // ── Taxonomy table ────────────────────────────────────────────────────
        Message::ShowTaxonomy   => { app.show_taxonomy = true;  Task::none() }
        Message::CloseTaxonomy  => { app.show_taxonomy = false; Task::none() }
        Message::SortByTaxonomy(level) => {
            if let Some(arc) = &mut app.document {
                let mut aln = (**arc).clone();
                sort_by_taxonomy(&mut aln, level);
                *arc = Arc::new(aln);
                app.status = if level == usize::MAX {
                    "Sorted by organism name".to_string()
                } else {
                    format!("Sorted by taxonomy rank {}", level + 1)
                };
            }
            Task::none()
        }

        Message::ExportSvgDialog => {
            if app.document.is_some() {
                Task::perform(save_svg_dialog(), Message::ExportSvgChosen)
            } else {
                Task::none()
            }
        }

        Message::ExportSvgChosen(None) => Task::none(),

        Message::ExportSvgChosen(Some(path)) => {
            if let Some(aln) = &app.document {
                let aln        = Arc::clone(aln);
                let rows       = app.view.selected.clone();
                let col_start  = app.view.selected_cols.iter().next().copied()
                    .unwrap_or(app.view.scroll_col);
                let col_end    = app.view.selected_cols.iter().last().copied()
                    .map(|c| c + 1)
                    .unwrap_or_else(|| aln.col_count());
                let color_table = Arc::clone(&app.color_table);
                app.status = format!("Exporting SVG…");
                Task::perform(
                    async move {
                        let svg = crate::export_svg::alignment_to_svg(
                            &aln, &rows, col_start, col_end, &color_table,
                        );
                        std::fs::write(&path, svg.as_bytes())
                            .map_err(|e| e.to_string())
                    },
                    Message::ExportSvgDone,
                )
            } else {
                Task::none()
            }
        }

        Message::ExportSvgDone(Ok(())) => {
            app.status = "SVG exported.".to_string();
            Task::none()
        }
        Message::ExportSvgDone(Err(e)) => {
            app.status = format!("SVG export failed: {e}");
            Task::none()
        }

        // ── Plasmid map viewer ────────────────────────────────────────────────

        Message::ShowPlasmid => {
            let idx = app.view.selected.iter().next().copied().unwrap_or(0);
            app.plasmid_seq_idx = idx;
            app.show_plasmid = true;
            app.plasmid_re_cache = None;
            if app.plasmid_show_re {
                if let Some(aln) = app.document.clone() {
                    return Task::perform(
                        compute_plasmid_re_sites(aln, idx),
                        Message::PlasmidReSitesDone,
                    );
                }
            }
            Task::none()
        }
        Message::ClosePlasmid => { app.show_plasmid = false; Task::none() }
        Message::PlasmidPrevSeq => {
            if let Some(aln) = &app.document {
                if app.plasmid_seq_idx > 0 {
                    app.plasmid_seq_idx -= 1;
                } else {
                    app.plasmid_seq_idx = aln.seq_count().saturating_sub(1);
                }
            }
            app.plasmid_selected_feat = None;
            app.plasmid_re_cache = None;
            if app.plasmid_show_re {
                if let Some(aln) = app.document.clone() {
                    let idx = app.plasmid_seq_idx;
                    return Task::perform(
                        compute_plasmid_re_sites(aln, idx),
                        Message::PlasmidReSitesDone,
                    );
                }
            }
            Task::none()
        }
        Message::PlasmidNextSeq => {
            if let Some(aln) = &app.document {
                let n = aln.seq_count();
                if n > 0 {
                    app.plasmid_seq_idx = (app.plasmid_seq_idx + 1) % n;
                }
            }
            app.plasmid_selected_feat = None;
            app.plasmid_re_cache = None;
            if app.plasmid_show_re {
                if let Some(aln) = app.document.clone() {
                    let idx = app.plasmid_seq_idx;
                    return Task::perform(
                        compute_plasmid_re_sites(aln, idx),
                        Message::PlasmidReSitesDone,
                    );
                }
            }
            Task::none()
        }
        Message::PlasmidToggleFeatures => {
            app.plasmid_show_features = !app.plasmid_show_features;
            Task::none()
        }
        Message::PlasmidToggleRe => {
            app.plasmid_show_re = !app.plasmid_show_re;
            if app.plasmid_show_re {
                app.plasmid_re_cache = None;
                if let Some(aln) = app.document.clone() {
                    let idx = app.plasmid_seq_idx;
                    return Task::perform(
                        compute_plasmid_re_sites(aln, idx),
                        Message::PlasmidReSitesDone,
                    );
                }
            }
            Task::none()
        }
        Message::PlasmidReSitesDone(sites) => {
            app.plasmid_re_cache = Some(Arc::new(sites));
            Task::none()
        }
        Message::ExportPlasmidSvg => {
            if app.document.is_some() {
                Task::perform(save_plasmid_svg_dialog(), Message::ExportPlasmidSvgChosen)
            } else {
                Task::none()
            }
        }
        Message::ExportPlasmidSvgChosen(None) => Task::none(),
        Message::ExportPlasmidSvgChosen(Some(path)) => {
            if let Some(aln) = &app.document {
                let aln            = Arc::clone(aln);
                let seq_idx        = app.plasmid_seq_idx;
                let show_features  = app.plasmid_show_features;
                let show_re        = app.plasmid_show_re;
                app.status = "Exporting plasmid SVG…".to_string();
                Task::perform(
                    async move {
                        let svg = crate::views::plasmid_to_svg(&aln, seq_idx, show_features, show_re);
                        std::fs::write(&path, svg.as_bytes()).map_err(|e| e.to_string())
                    },
                    Message::ExportPlasmidSvgDone,
                )
            } else {
                Task::none()
            }
        }
        Message::ExportPlasmidSvgDone(Ok(())) => {
            app.status = "Plasmid SVG exported.".to_string();
            Task::none()
        }
        Message::ExportPlasmidSvgDone(Err(e)) => {
            app.status = format!("Plasmid SVG export failed: {e}");
            Task::none()
        }

        Message::PlasmidSelectFeature(idx) => {
            app.plasmid_selected_feat = idx;
            // Populate edit buffers from the selected feature.
            if let (Some(i), Some(aln)) = (idx, &app.document) {
                if let Some(seq) = aln.sequences.get(app.plasmid_seq_idx) {
                    if let Some(feat) = seq.features.get(i) {
                        app.plasmid_edit_name  = feat.name.clone();
                        app.plasmid_edit_start = (feat.start + 1).to_string();
                        app.plasmid_edit_end   = (feat.end   + 1).to_string();
                        app.plasmid_edit_color = feat.color;
                        app.plasmid_edit_hex   = format!("#{:02X}{:02X}{:02X}",
                            (feat.color.r * 255.0).round() as u8,
                            (feat.color.g * 255.0).round() as u8,
                            (feat.color.b * 255.0).round() as u8);
                    }
                }
            }
            Task::none()
        }

        Message::PlasmidEditName(s)  => { app.plasmid_edit_name  = s; Task::none() }
        Message::PlasmidEditStart(s) => { app.plasmid_edit_start = s; Task::none() }
        Message::PlasmidEditEnd(s)   => { app.plasmid_edit_end   = s; Task::none() }
        Message::PlasmidEditColor(c) => {
            app.plasmid_edit_color = c;
            app.plasmid_edit_hex = format!("#{:02X}{:02X}{:02X}",
                (c.r * 255.0).round() as u8,
                (c.g * 255.0).round() as u8,
                (c.b * 255.0).round() as u8);
            Task::none()
        }
        Message::PlasmidEditHex(s) => {
            // Update the color if the hex parses as a valid #RRGGBB; keep the string regardless.
            if let Some(c) = parse_hex_color(&s) {
                app.plasmid_edit_color = c;
            }
            app.plasmid_edit_hex = s;
            Task::none()
        }

        Message::PlasmidApplyEdit => {
            if let (Some(idx), Some(arc)) = (app.plasmid_selected_feat, &mut app.document) {
                let seq_idx   = app.plasmid_seq_idx;
                let new_start = app.plasmid_edit_start.parse::<usize>()
                    .map(|v| v.saturating_sub(1))
                    .unwrap_or_else(|_| {
                        (**arc).sequences.get(seq_idx)
                            .and_then(|s| s.features.get(idx))
                            .map(|f| f.start)
                            .unwrap_or(0)
                    });
                let new_end = app.plasmid_edit_end.parse::<usize>()
                    .map(|v| v.saturating_sub(1))
                    .unwrap_or_else(|_| {
                        (**arc).sequences.get(seq_idx)
                            .and_then(|s| s.features.get(idx))
                            .map(|f| f.end)
                            .unwrap_or(0)
                    });
                let cmd = Box::new(EditFeature::new(
                    seq_idx, idx,
                    app.plasmid_edit_name.clone(),
                    new_start, new_end,
                    app.plasmid_edit_color,
                ));
                let mut aln = (**arc).clone();
                app.history.execute(cmd, &mut aln);
                *arc = Arc::new(aln);
            }
            Task::none()
        }

        Message::PlasmidDeleteFeature(idx) => {
            if let Some(arc) = &mut app.document {
                let seq_idx = app.plasmid_seq_idx;
                let cmd = Box::new(DeleteFeatureCmd::new(seq_idx, idx));
                let mut aln = (**arc).clone();
                app.history.execute(cmd, &mut aln);
                *arc = Arc::new(aln);
                app.plasmid_selected_feat = None;
            }
            Task::none()
        }

        Message::PlasmidAddFeature => {
            use helixview_core::feature::Feature;
            if let Some(arc) = &mut app.document {
                let seq_idx = app.plasmid_seq_idx;
                let (new_idx, end_pos) = if let Some(seq) = (**arc).sequences.get(seq_idx) {
                    (seq.features.len(), seq.residues.len().min(100))
                } else {
                    (0, 100)
                };
                let mut feat = Feature::new("New feature", 0, end_pos.saturating_sub(1));
                // Cycle through distinct colors so successive features are visually distinct.
                const NEW_FEAT_PALETTE: [(u8,u8,u8); 8] = [
                    ( 52, 120, 205), (220, 115,  25), ( 50, 180,  65),
                    (215,  50,  50), (140,  50, 190), (  5, 165, 165),
                    (215,  50, 120), (189, 166,  12),
                ];
                let (r, g, b) = NEW_FEAT_PALETTE[new_idx % NEW_FEAT_PALETTE.len()];
                feat.color = helixview_core::color::Color::from_u8(r, g, b);
                let cmd  = Box::new(AddFeature { seq_idx, feature: feat });
                let mut aln = (**arc).clone();
                app.history.execute(cmd, &mut aln);
                *arc = Arc::new(aln);
                app.plasmid_selected_feat = Some(new_idx);
                app.plasmid_edit_name  = "New feature".to_string();
                app.plasmid_edit_start = "1".to_string();
                app.plasmid_edit_end   = end_pos.to_string();
            }
            Task::none()
        }

        Message::PlasmidZoomDelta(d) => {
            app.plasmid_zoom = (app.plasmid_zoom * (1.0 + d)).clamp(0.2, 5.0);
            Task::none()
        }
        Message::PlasmidZoomIn => {
            app.plasmid_zoom = (app.plasmid_zoom + 0.15).min(5.0);
            Task::none()
        }
        Message::PlasmidZoomOut => {
            app.plasmid_zoom = (app.plasmid_zoom - 0.15).max(0.2);
            Task::none()
        }
        Message::PlasmidZoomReset => {
            app.plasmid_zoom  = 1.0;
            app.plasmid_pan_x = 0.0;
            app.plasmid_pan_y = 0.0;
            Task::none()
        }
        Message::PlasmidPanDelta(dx, dy) => {
            app.plasmid_pan_x += dx;
            app.plasmid_pan_y += dy;
            Task::none()
        }

        // ── Tab management ────────────────────────────────────────────────────

        Message::NewTab => {
            app.open_in_tab(SeqAlignment::new("Untitled"));
            app.status = "New tab opened".to_string();
            Task::none()
        }

        Message::SwitchTab(idx) => {
            app.switch_to(idx);
            Task::none()
        }

        Message::CloseTab(idx) => {
            if app.tabs.len() <= 1 {
                // Only one tab — just clear it rather than removing.
                app.document    = None;
                app.view        = ViewState::default();
                app.history     = EditHistory::new(100);
                app.show_restr_map = false;
                app.protein_view_doc = None;
                app.tabs[0]     = TabRecord::empty();
                app.active_tab  = 0;
                app.status      = "Closed tab".to_string();
            } else {
                // Switch away first if we're closing the active tab.
                if idx == app.active_tab {
                    let go_to = if idx + 1 < app.tabs.len() { idx + 1 } else { idx - 1 };
                    app.switch_to(go_to);
                }
                // Adjust active_tab index after removal.
                let remove_at = if idx < app.active_tab { idx } else { idx };
                app.tabs.remove(remove_at);
                if app.active_tab > 0 && remove_at < app.active_tab {
                    app.active_tab -= 1;
                }
                app.status = "Closed tab".to_string();
            }
            Task::none()
        }

        Message::TreeReorderAlignment => {
            if let (Some(tree), Some(doc)) = (&app.tree_data, &mut app.document) {
                let names: Vec<String> = doc.sequences.iter().map(|s| s.name.clone()).collect();
                let leaves = tree.leaf_names();
                let matched = leaves.iter()
                    .filter(|leaf| names.iter().any(|n| {
                        let lu = leaf.to_ascii_uppercase();
                        let nu = n.to_ascii_uppercase();
                        nu == lu || nu.starts_with(&lu) || lu.starts_with(&nu)
                    }))
                    .count();
                let order = tree.alignment_order(&names);
                let mut aln = (**doc).clone();
                let old_seqs = aln.sequences.clone();
                for (new_pos, &old_pos) in order.iter().enumerate() {
                    aln.sequences[new_pos] = old_seqs[old_pos].clone();
                }
                *doc = Arc::new(aln);
                if matched == 0 {
                    app.status = format!(
                        "Reorder: 0/{} leaves matched — check that tree taxon names match sequence names.",
                        leaves.len()
                    );
                } else {
                    app.status = format!(
                        "Reordered: {matched}/{} leaves matched. Click < Back to see alignment.",
                        leaves.len()
                    );
                }
            } else if app.document.is_none() {
                app.status = "No alignment open — open a file first, then reorder.".to_string();
            }
            Task::none()
        }

        Message::AppCloseRequested(id) => {
            auto_save_session(app);
            iced::window::close(id)
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Build a protein view alignment from a nucleotide alignment using the given frame (0–5).
fn make_protein_alignment(aln: &SeqAlignment, frame: usize) -> SeqAlignment {
    let frame_label: i8 = if frame < 3 { (frame as i8) + 1 } else { -((frame as i8) - 2) };
    let mut out = SeqAlignment::new(&format!("{} [Protein {:+}]", aln.name, frame_label));
    for seq in &aln.sequences {
        if !seq.seq_type.is_nucleic() { continue; }
        let protein = helixview_analysis::translate_frame(&seq.residues, frame);
        let mut ps = seq.clone();
        ps.residues  = protein;
        ps.gap_locks.clear();
        ps.seq_type  = helixview_core::SequenceType::Protein;
        out.push(ps);
    }
    // Pad to equal width.
    let max_len = out.sequences.iter().map(|s| s.residues.len()).max().unwrap_or(0);
    for seq in &mut out.sequences {
        seq.residues.resize(max_len, b'-');
        seq.gap_locks.resize(max_len, false);
    }
    out
}

// ── Subscription ─────────────────────────────────────────────────────────────

pub fn subscription(_app: &HelixViewApp) -> iced::Subscription<Message> {
    use iced::event;
    use iced::keyboard::Event as KbEvent;

    iced::Subscription::batch([
        keyboard::on_key_press(|key, mods| {
            match key {
                Key::Named(Named::ArrowLeft)  => Some(Message::GridScrolled { dx: -1, dy:   0 }),
                Key::Named(Named::ArrowRight) => Some(Message::GridScrolled { dx:  1, dy:   0 }),
                Key::Named(Named::ArrowUp)    => Some(Message::GridScrolled { dx:  0, dy:  -1 }),
                Key::Named(Named::ArrowDown)  => Some(Message::GridScrolled { dx:  0, dy:   1 }),
                Key::Named(Named::PageUp)     => Some(Message::GridScrolled { dx:  0, dy: -20 }),
                Key::Named(Named::PageDown)   => Some(Message::GridScrolled { dx:  0, dy:  20 }),
                Key::Named(Named::Home)       => Some(Message::ScrollColSet(0)),
                Key::Named(Named::End)        => Some(Message::ScrollToEnd),
                // Backspace/Delete/Space routing is handled by update() based on edit mode.
                Key::Named(Named::Backspace)  => Some(Message::BackspacePressed { shift: mods.shift() }),
                Key::Named(Named::Delete)     => Some(Message::DeleteSelectedColumns),
                Key::Named(Named::Space)      => Some(Message::SpacePressed { shift: mods.shift() }),
                Key::Named(Named::Tab)        => Some(Message::GridScrolled { dx: 1, dy: 0 }),
                _ => {
                    if mods.control() {
                        match key {
                            Key::Character(c) if c.as_str() == "z" => Some(Message::Undo),
                            Key::Character(c) if c.as_str() == "y" => Some(Message::Redo),
                            Key::Character(c) if c.as_str() == "s" => Some(Message::SaveFileDialog),
                            Key::Character(c) if c.as_str() == "o" => Some(Message::OpenFileDialog),
                            Key::Character(c) if c.as_str() == "n" => Some(Message::NewAlignment),
                            Key::Character(c) if c.as_str() == "f" => Some(Message::ToggleSearch),
                            Key::Character(c) if c.as_str() == "a" => Some(Message::SelectAll),
                            Key::Character(c) if c.as_str() == "d" => Some(Message::SelectNone),
                            Key::Character(c) if c.as_str() == "=" || c.as_str() == "+" => Some(Message::ZoomIn),
                            Key::Character(c) if c.as_str() == "-" => Some(Message::ZoomOut),
                            Key::Character(c) if c.as_str() == "0" => Some(Message::ZoomReset),
                            Key::Character(c) if c.as_str() == "p" => Some(Message::OpenCommandPalette),
                            _ => None,
                        }
                    } else {
                        match key {
                            Key::Named(Named::Escape) => Some(Message::SearchClose),
                            Key::Named(Named::ArrowUp)   => Some(Message::PaletteUp),
                            Key::Named(Named::ArrowDown) => Some(Message::PaletteDown),
                            Key::Character(c) => {
                                // Single printable char → residue editing if a cell is selected.
                                // (The app update handler will ignore this when no cell is selected.)
                                let s = c.as_str();
                                if s.len() == 1 {
                                    let ch = s.chars().next().unwrap();
                                    if ch.is_ascii_alphanumeric() || ch == '-' || ch == '.' {
                                        Some(Message::ResidueTyped(ch))
                                    } else {
                                        None
                                    }
                                } else {
                                    None
                                }
                            }
                            _ => None,
                        }
                    }
                }
            }
        }),
        event::listen_with(|ev, _status, _window| {
            if let iced::Event::Keyboard(KbEvent::ModifiersChanged(mods)) = ev {
                Some(Message::ModifiersChanged(mods))
            } else {
                None
            }
        }),
    ])
}

// ── View ──────────────────────────────────────────────────────────────────────

fn tab_bar(app: &HelixViewApp) -> iced::widget::Row<'_, Message> {
    use iced::widget::{button, row, text, horizontal_space};
    use iced::Length;

    let mut r = row![].spacing(2).padding([2, 4]);
    for (i, rec) in app.tabs.iter().enumerate() {
        let active = i == app.active_tab;
        // Active tab's live document is in app.document, not the stale snapshot.
        let name = if active {
            app.document.as_ref()
                .map(|d| d.name.clone())
                .filter(|n| !n.is_empty())
                .unwrap_or_else(|| "New Tab".to_string())
        } else {
            HelixViewApp::tab_name(rec)
        };
        let tab_btn = button(text(name).size(11))
            .padding([3, 10])
            .style(if active { button::primary } else { button::secondary })
            .on_press(Message::SwitchTab(i));
        let close_btn = button(text("×").size(11))
            .padding([3, 5])
            .style(button::secondary)
            .on_press(Message::CloseTab(i));
        r = r.push(tab_btn).push(close_btn);
    }
    r = r.push(horizontal_space());
    r = r.push(
        button(text("+").size(11))
            .padding([3, 8])
            .style(button::secondary)
            .on_press(Message::NewTab),
    );
    r
}

pub fn view(app: &HelixViewApp) -> Element<'_, Message> {
    use iced::widget::column;
    use iced::widget::container;
    use iced::Length;

    let content = view_content(app);
    let tabs = tab_bar(app);

    // Thin separator line below tab bar.
    let sep = container(iced::widget::horizontal_space())
        .width(Length::Fill)
        .height(iced::Length::Fixed(1.0))
        .style(|_| iced::widget::container::Style {
            background: Some(iced::Background::Color(iced::Color::from_rgb(0.75, 0.76, 0.80))),
            ..Default::default()
        });

    column![tabs, sep, content]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

fn view_content(app: &HelixViewApp) -> Element<'_, Message> {
    // Color table editor.
    if app.show_color_editor {
        return color_editor_view(&app.color_table, app.color_editor_selected);
    }
    // Conservation search view.
    if app.show_conservation {
        if let Some(aln) = &app.document {
            return conservation_view(aln, app.conservation_threshold, app.conservation_min_width);
        }
    }
    // Sequence raw editor takes top priority.
    if let (Some(idx), Some(content)) = (app.seq_editor_idx, &app.seq_editor_content) {
        let seq_name = app.document.as_ref()
            .and_then(|aln| aln.sequences.get(idx))
            .map(|s| s.name.as_str())
            .unwrap_or("?");
        return seq_editor_view(seq_name, content, idx);
    }
    // ORF finder results take priority.
    if let Some((name, orfs)) = &app.orf_results {
        return crate::views::orfs_view(name, orfs);
    }
    // Six-frame translation view.
    if let Some((name, nt_len, frames)) = &app.six_frame_result {
        return crate::views::sixframe_view(name, *nt_len, frames);
    }
    // Pairwise alignment results.
    if let Some(result) = &app.pairwise_result {
        return crate::views::pairwise_view(result);
    }
    // Identity matrix view.
    if app.show_identity_matrix {
        if let Some(aln) = &app.document {
            return identity_matrix_view(aln, app.identity_matrix_data.as_deref());
        }
    }
    // Positional column summary.
    if app.show_col_summary {
        if let Some(aln) = &app.document {
            return col_summary_view(aln);
        }
    }
    // Dot-plot view.
    if app.show_dot_plot {
        if let Some(aln) = &app.document {
            return dot_plot_view(
                aln,
                app.dot_plot_seq_a,
                app.dot_plot_seq_b,
                app.dot_plot_window,
                app.dot_plot_min_match,
            );
        }
    }
    // Tree viewer.
    if let Some(tree) = &app.tree_data {
        return crate::views::tree_view(tree, &app.tree_selected);
    }
    // ABI trace view.
    if let Some(trace) = &app.trace_data {
        return crate::views::trace_view(trace, app.trace_scroll, app.trace_zoom, app.trace_quant_ch_a, app.trace_quant_ch_b);
    }
    // Accessory app manager.
    if app.show_accessories {
        return crate::views::accessories_view(
            &app.accessories,
            app.editing_accessory,
            &app.acc_edit_name,
            &app.acc_edit_cmd,
            app.document.is_some(),
        );
    }
    // Mutual information view.
    if app.show_mutual_info {
        if let Some(aln) = &app.document {
            return crate::views::mutual_info_view(
                app.mi_result.as_deref(),
                app.mi_min_obs,
                app.mi_top_n,
                app.mi_stem_min_wc_types,
                app.mi_stem_min_wc_frac,
                &app.mi_tab,
                app.mi_running,
            );
        }
    }
    // Composition analysis view.
    if app.show_composition {
        if let Some(aln) = &app.document {
            return crate::views::composition_view(aln);
        }
    }
    // Hydrophobicity profile view.
    if app.show_hydrophobicity {
        if let Some(aln) = &app.document {
            return crate::views::hydrophobicity_view(aln, app.hydro_window);
        }
    }
    // Tm calculator view.
    if app.show_oligo_tm {
        if let Some(aln) = &app.document {
            return crate::views::oligo_tm_view(aln);
        }
    }
    // Help / shortcuts view.
    if app.show_help {
        return crate::views::help_view();
    }
    // Preferences dialog.
    if app.show_prefs {
        return crate::views::prefs_dialog(&app.prefs);
    }
    // Shaded export dialog.
    if app.show_shaded_export {
        return crate::views::shaded_export_dialog(&app.shaded_export_opts);
    }
    // Taxonomy table view.
    if app.show_taxonomy {
        if let Some(aln) = &app.document {
            return crate::views::taxonomy_view(aln);
        }
    }
    // Plasmid map view.
    if app.show_plasmid {
        if let Some(aln) = &app.document {
            let re_cache = app.plasmid_re_cache.clone();
            return crate::views::plasmid_view(
                aln,
                app.plasmid_seq_idx,
                app.plasmid_show_features,
                app.plasmid_show_re,
                app.plasmid_selected_feat,
                &app.plasmid_edit_name,
                &app.plasmid_edit_start,
                &app.plasmid_edit_end,
                app.plasmid_edit_color,
                &app.plasmid_edit_hex,
                app.plasmid_zoom,
                app.plasmid_pan_x,
                app.plasmid_pan_y,
                re_cache,
            );
        }
    }
    // Text alignment export view.
    if app.show_text_export {
        if let Some(aln) = &app.document {
            let preview = crate::text_export::format_alignment(aln, &app.text_export_opts);
            return crate::views::text_export_dialog(&app.text_export_opts, &preview);
        }
    }
    // BLAST results view.
    if !matches!(app.blast, BlastState::Idle) {
        return crate::views::blast_view(&app.blast, None);
    }
    let base: Element<Message> = match &app.document {
        None => welcome_view(),
        Some(aln) => {
            // When protein view is active, display the derived protein alignment.
            let display = app.protein_view_doc.as_ref().unwrap_or(aln);
            let pairing_pairs: Vec<crate::widgets::PairingPair> = if app.show_pairing_arcs {
                app.mi_result.as_deref()
                    .map(|mi| mi.top_covarying.iter()
                        .filter(|p| p.cov_score >= app.pairing_arc_threshold)
                        .take(60)
                        .map(|p| (p.col_a, p.col_b, p.cov_score))
                        .collect())
                    .unwrap_or_default()
            } else { vec![] };
            let pairing_max = app.mi_result.as_deref()
                .map(|mi| mi.top_covarying.first().map_or(1.0, |p| p.cov_score.max(1.0)))
                .unwrap_or(1.0);
            let mi_available = app.mi_result.is_some();
            let mi_stale = mi_available && app.mi_aln_ptr
                .map_or(false, |p| app.document.as_ref().map_or(true, |a| Arc::as_ptr(a) as usize != p));
            alignment_view(display, &app.view, &app.status, app.title_menu, &app.rename_input, app.show_restr_map, &app.color_table, app.computing.as_deref(), pairing_pairs, pairing_max, app.show_pairing_arcs, app.pairing_arc_threshold, mi_available, mi_stale)
        }
    };
    // Command palette overlays the base view.
    if app.show_command_palette {
        use iced::widget::stack;
        return stack![
            base,
            crate::views::command_palette(&app.palette_query, app.palette_selected),
        ].into();
    }
    base
}

// ── Async helpers ──────────────────────────────────────────────────────────────

fn parse_hex_color(s: &str) -> Option<helixview_core::color::Color> {
    let s = s.trim().trim_start_matches('#');
    if s.len() != 6 { return None; }
    let r = u8::from_str_radix(&s[0..2], 16).ok()?;
    let g = u8::from_str_radix(&s[2..4], 16).ok()?;
    let b = u8::from_str_radix(&s[4..6], 16).ok()?;
    Some(helixview_core::color::Color::from_u8(r, g, b))
}

async fn pick_file() -> Option<std::path::PathBuf> {
    rfd::AsyncFileDialog::new()
        .set_title("Open sequence file")
        .add_filter("Sequence files",
            &["fasta", "fa", "fna", "faa", "ffn", "frn",
              "gb", "gbk", "genbank", "bio",
              "aln", "pir", "nbrf",
              "phy", "phylip", "nex", "nexus", "nxs", "msf",
              "embl", "dat",
              "sto", "stockholm", "stk"])
        .add_filter("All files", &["*"])
        .pick_file()
        .await
        .map(|h| h.path().to_path_buf())
}

async fn save_file_dialog() -> Option<std::path::PathBuf> {
    rfd::AsyncFileDialog::new()
        .set_title("Save alignment as…")
        .add_filter("FASTA",        &["fasta", "fa", "fna", "faa"])
        .add_filter("ClustalW",     &["aln"])
        .add_filter("NBRF/PIR",     &["pir", "nbrf"])
        .add_filter("Phylip",       &["phy", "phylip"])
        .add_filter("NEXUS",        &["nex", "nexus"])
        .add_filter("MSF (GCG)",    &["msf"])
        .add_filter("EMBL",         &["embl"])
        .add_filter("All files",    &["*"])
        .save_file()
        .await
        .map(|h| h.path().to_path_buf())
}

async fn save_svg_dialog() -> Option<std::path::PathBuf> {
    rfd::AsyncFileDialog::new()
        .set_title("Export alignment as SVG…")
        .add_filter("SVG", &["svg"])
        .add_filter("All files", &["*"])
        .save_file()
        .await
        .map(|h| h.path().to_path_buf())
}

async fn pick_tree_file() -> Option<std::path::PathBuf> {
    rfd::AsyncFileDialog::new()
        .set_title("Open phylogenetic tree file")
        .add_filter("Tree files", &["nwk", "newick", "nex", "nexus", "nxs", "tree", "tre"])
        .add_filter("All files", &["*"])
        .pick_file()
        .await
        .map(|h| h.path().to_path_buf())
}

async fn pick_trace_file() -> Option<std::path::PathBuf> {
    rfd::AsyncFileDialog::new()
        .set_title("Open ABI trace file")
        .add_filter("ABI trace", &["ab1", "abi", "AB1", "ABI"])
        .add_filter("All files", &["*"])
        .pick_file()
        .await
        .map(|h| h.path().to_path_buf())
}

async fn save_project_dialog() -> Option<std::path::PathBuf> {
    rfd::AsyncFileDialog::new()
        .set_title("Save HelixView project…")
        .add_filter("HelixView project", &["hvproj"])
        .add_filter("All files", &["*"])
        .save_file()
        .await
        .map(|h| h.path().to_path_buf())
}

async fn open_project_dialog() -> Option<std::path::PathBuf> {
    rfd::AsyncFileDialog::new()
        .set_title("Open HelixView project…")
        .add_filter("HelixView project", &["hvproj"])
        .add_filter("All files", &["*"])
        .pick_file()
        .await
        .map(|h| h.path().to_path_buf())
}

async fn save_shaded_svg_dialog() -> Option<std::path::PathBuf> {
    rfd::AsyncFileDialog::new()
        .set_title("Export shaded alignment graphic…")
        .add_filter("SVG", &["svg"])
        .add_filter("All files", &["*"])
        .save_file()
        .await
        .map(|h| h.path().to_path_buf())
}

async fn save_plasmid_svg_dialog() -> Option<std::path::PathBuf> {
    rfd::AsyncFileDialog::new()
        .set_title("Export plasmid map as SVG…")
        .add_filter("SVG", &["svg"])
        .add_filter("All files", &["*"])
        .save_file()
        .await
        .map(|h| h.path().to_path_buf())
}

async fn save_text_dialog() -> Option<std::path::PathBuf> {
    rfd::AsyncFileDialog::new()
        .set_title("Export alignment as text…")
        .add_filter("Text", &["txt"])
        .add_filter("All files", &["*"])
        .save_file()
        .await
        .map(|h| h.path().to_path_buf())
}

async fn save_pdf_dialog() -> Option<std::path::PathBuf> {
    rfd::AsyncFileDialog::new()
        .set_title("Export alignment as PDF…")
        .add_filter("PDF", &["pdf"])
        .add_filter("All files", &["*"])
        .save_file()
        .await
        .map(|h| h.path().to_path_buf())
}

async fn export_pdf_task(
    aln:         Arc<SeqAlignment>,
    opts:        ShadedExportOptions,
    col_start:   usize,
    col_end:     usize,
    color_table: Arc<ColorTable>,
    path:        std::path::PathBuf,
) -> Result<(), String> {
    tokio::task::spawn_blocking(move || {
        let all_rows: std::collections::BTreeSet<usize> =
            (0..aln.seq_count()).collect();
        let svg = crate::export_shaded::shaded_alignment_svg(
            &aln,
            &all_rows,
            col_start,
            col_end,
            &color_table,
            &opts.scheme,
            opts.threshold,
            opts.show_ruler,
            opts.show_consensus,
        );
        let mut usvg_opts = svg2pdf::usvg::Options::default();
        usvg_opts.fontdb_mut().load_system_fonts();
        let tree = svg2pdf::usvg::Tree::from_str(&svg, &usvg_opts)
            .map_err(|e| e.to_string())?;
        let pdf = svg2pdf::to_pdf(
            &tree,
            svg2pdf::ConversionOptions::default(),
            svg2pdf::PageOptions::default(),
        )
        .map_err(|e| e.to_string())?;
        std::fs::write(&path, pdf).map_err(|e| e.to_string())
    })
    .await
    .map_err(|e| e.to_string())?
}

async fn load_file(path: std::path::PathBuf) -> Result<SeqAlignment, String> {
    tokio::task::spawn_blocking(move || {
        helixview_formats::open(&path).map_err(|e: helixview_formats::FormatError| e.to_string())
    })
    .await
    .map_err(|e: tokio::task::JoinError| e.to_string())?
}

/// Compute restriction-enzyme sites for the plasmid sequence on a background thread.
/// Returns (enzyme_name, 1-based_position) pairs sorted by position.
async fn compute_plasmid_re_sites(
    aln: Arc<SeqAlignment>,
    seq_idx: usize,
) -> Vec<helixview_analysis::RestrictionSite> {
    tokio::task::spawn_blocking(move || {
        let seq = match aln.sequences.get(seq_idx) {
            Some(s) => s,
            None    => return vec![],
        };
        let raw: Vec<u8> = seq.residues.iter()
            .filter(|&&b| !matches!(b, b'-' | b'.' | b'~'))
            .map(|&b| b.to_ascii_uppercase())
            .collect();
        let mut sites = helixview_analysis::find_sites(&raw);
        sites.sort_by_key(|s| s.true_pos);
        sites
    })
    .await
    .unwrap_or_default()
}

// ── Entry point ───────────────────────────────────────────────────────────────

/// Build a 32×32 RGBA icon: dark-green background, white "H" letterform.
fn make_icon() -> Option<iced::window::Icon> {
    const S: usize = 32;
    let mut px = vec![0u8; S * S * 4];

    // Background: brand green (26, 133, 71)
    for i in 0..S * S {
        px[i * 4]     = 26;
        px[i * 4 + 1] = 133;
        px[i * 4 + 2] = 71;
        px[i * 4 + 3] = 255;
    }

    // Draw a white "H": two vertical bars (cols 7-11, 20-24) + crossbar (rows 13-18)
    let paint = |px: &mut Vec<u8>, col: usize, row: usize| {
        if col < S && row < S {
            let i = (row * S + col) * 4;
            px[i] = 255; px[i+1] = 255; px[i+2] = 255; px[i+3] = 255;
        }
    };

    for row in 5..27usize {
        for col in 7..12usize  { paint(&mut px, col, row); } // left bar
        for col in 20..25usize { paint(&mut px, col, row); } // right bar
    }
    for row in 13..19usize {
        for col in 12..20usize { paint(&mut px, col, row); } // crossbar
    }

    iced::window::icon::from_rgba(px, S as u32, S as u32).ok()
}

pub fn run() -> iced::Result {
    tracing_subscriber::fmt::init();

    let window = {
        let mut w = iced::window::Settings::default();
        w.icon = make_icon();
        // Set the Wayland app-id / X11 WM_CLASS so the compositor shows
        // "helixview" instead of the winit default (which renders as "W").
        w.platform_specific.application_id = "helixview".to_string();
        w
    };

    iced::application(app_title, update, view)
        .subscription(subscription)
        .theme(|_| Theme::Light)
        .window_size((1200.0, 780.0))
        .antialiasing(false)
        .window(window)
        .run_with(|| {
            let app = HelixViewApp::default();
            let task = if let Some(path) = crate::prefs::session_path() {
                if path.exists() {
                    Task::done(Message::OpenProjectChosen(Some(path)))
                } else {
                    Task::none()
                }
            } else {
                Task::none()
            };
            (app, task)
        })
}

// ── Utility ───────────────────────────────────────────────────────────────────

/// Silently save the current session to the auto-session path.
/// Skips if the alignment is large (would block the main thread for too long).
fn auto_save_session(app: &HelixViewApp) {
    // Skip auto-save for large alignments — serialising millions of residues
    // to JSON on the main thread causes multi-second UI freezes.
    let total_residues: usize = app.document
        .as_ref()
        .map_or(0, |a| a.seq_count().saturating_mul(a.col_count()));
    if total_residues > 500_000 {
        return;
    }

    if let Some(path) = crate::prefs::session_path() {
        // Serialise on the main thread (small alignment, fast), then write async.
        if let Ok(json) = crate::project::to_json(app) {
            tokio::task::spawn(async move {
                let _ = tokio::fs::write(&path, json).await;
            });
        }
    }
}

fn app_title(app: &HelixViewApp) -> String {
    match &app.document {
        None => "HelixView".to_string(),
        Some(aln) if !aln.name.is_empty() => format!("{} — HelixView", aln.name),
        _ => "HelixView".to_string(),
    }
}

/// Returns true if residue `r` matches IUPAC ambiguity code `p`.
/// Both bytes should be uppercase ASCII. `p` = pattern byte, `r` = residue byte.
fn iupac_matches(p: u8, r: u8) -> bool {
    if p == r { return true; }
    // IUPAC nucleotide ambiguity codes and wildcards.
    // For amino acid letters that aren't ambiguity codes (F, E, I, L, P, Q, W, Y etc.)
    // the p == r check above handles exact matching; no special expansion needed.
    match p {
        b'N' => matches!(r, b'A'|b'C'|b'G'|b'T'|b'U'|b'N'),
        // X matches any non-gap residue (used in both nucleotide and protein searches)
        b'X' => r.is_ascii_alphabetic(),
        b'R' => matches!(r, b'A'|b'G'),
        b'Y' => matches!(r, b'C'|b'T'|b'U'),
        b'S' => matches!(r, b'G'|b'C'),
        b'W' => matches!(r, b'A'|b'T'|b'U'),
        b'K' => matches!(r, b'G'|b'T'|b'U'),
        b'M' => matches!(r, b'A'|b'C'),
        b'B' => matches!(r, b'C'|b'G'|b'T'|b'U'),
        b'D' => matches!(r, b'A'|b'G'|b'T'|b'U'),
        b'H' => matches!(r, b'A'|b'C'|b'T'|b'U'),
        b'V' => matches!(r, b'A'|b'C'|b'G'),
        // Any other letter (amino acid residues): exact match only (handled above).
        _ => false,
    }
}

/// Reorder sequences in `aln` by taxonomy rank or organism name.
/// `level == usize::MAX` → sort by organism name.
/// Otherwise sort by the `level`-th semicolon-delimited lineage rank.
fn sort_by_taxonomy(aln: &mut helixview_core::Alignment, level: usize) {
    fn parse_source_key(source: &str, level: usize) -> String {
        let mut in_lineage = false;
        let mut lineage_buf = String::new();
        let mut organism = String::new();
        for line in source.lines() {
            if let Some(rest) = line.strip_prefix("ORGANISM") {
                organism = rest.trim().to_string();
                in_lineage = true;
            } else if in_lineage {
                if !lineage_buf.is_empty() { lineage_buf.push(' '); }
                lineage_buf.push_str(line.trim());
            }
        }
        if level == usize::MAX {
            return if organism.is_empty() {
                source.lines().next().unwrap_or("").trim().to_string()
            } else {
                organism
            };
        }
        lineage_buf
            .split(';')
            .nth(level)
            .unwrap_or("")
            .trim()
            .trim_end_matches('.')
            .to_string()
    }

    aln.sequences.sort_by(|a, b| {
        let ka = a.genbank.source.as_deref()
            .map(|s| parse_source_key(s, level))
            .unwrap_or_default();
        let kb = b.genbank.source.as_deref()
            .map(|s| parse_source_key(s, level))
            .unwrap_or_default();
        ka.cmp(&kb)
    });
}

/// Add a signed delta to a scroll position, clamping to `0..len`.
fn clamp_add(pos: usize, delta: isize, len: usize) -> usize {
    if len == 0 { return 0; }
    ((pos as isize + delta).max(0) as usize).min(len - 1)
}
