/** @type {import('tailwindcss').Config} */

// Semantic names, so a component says what a colour is for rather than which
// shade it happens to be. The values come from CSS variables that `.dark` on
// <html> swaps — see src/app.css. Written as bare channels so Tailwind's
// opacity modifiers (bg-surface/80) still work.
const themed = name => `rgb(var(--c-${name}) / <alpha-value>)`;

export default {
  content: ['./src/**/*.{svelte,js,ts}', './index.html'],
  darkMode: 'class',
  theme: {
    extend: {
      fontFamily: {
        sans: ['Inter', 'system-ui', 'sans-serif'],
      },
      colors: {
        base: themed('base'),        // the page behind everything
        surface: themed('surface'),  // panels, bars, drawers
        sunken: themed('sunken'),    // recessed wells: code blocks, tracks
        raised: themed('raised'),    // controls sitting on a panel
        raised2: themed('raised-2'), // the same control, hovered or active
        strong: themed('strong'),    // headings
        fg: themed('fg'),            // body text
        fg2: themed('fg-2'),
        muted: themed('muted'),      // labels
        faint: themed('faint'),      // asides, units, placeholder text
        line: themed('line'),        // panel edges
        line2: themed('line-2'),     // an edge that has to be seen
        accent: themed('accent'),    // the active tab, the sort arrow
        link: themed('link'),
        link2: themed('link-2'),     // a link, hovered
        warn: themed('warn'),        // a figure the reader should not like
        caution: themed('caution'),
      },
    },
  },
  plugins: [],
};
