// Latin only. The unqualified imports pull every subset Inter ships — Greek,
// Cyrillic, Vietnamese and their extended forms — which is 56 font files and
// most of the bundle, to render taxonomy and sample identifiers that are Latin
// script by construction. latin-ext covers the diacritics a place name or an
// author's name can carry.
import '@fontsource/inter/latin-400.css';
import '@fontsource/inter/latin-500.css';
import '@fontsource/inter/latin-600.css';
import '@fontsource/inter/latin-700.css';
import '@fontsource/inter/latin-ext-400.css';
import '@fontsource/inter/latin-ext-500.css';
import '@fontsource/inter/latin-ext-600.css';
import '@fontsource/inter/latin-ext-700.css';
import './app.css';
import { mount } from 'svelte';
import App from './App.svelte';

const app = mount(App, { target: document.getElementById('app') });

export default app;
