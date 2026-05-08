import clsx from 'clsx';
import Link from '@docusaurus/Link';
import Layout from '@theme/Layout';
import useBaseUrl from '@docusaurus/useBaseUrl';
import styles from './index.module.css';

export default function Home(): JSX.Element {
  const logoUrl = useBaseUrl('/img/cbicall-logo.png');

  return (
    <Layout
      title="CBIcall"
      description="CNAG Biomedical Informatics framework for variant calling">
      <main>
        <section className={styles.hero}>
          <div className={styles.heroInner}>
            <img className={styles.heroLogo} src={logoUrl} alt="CBIcall" />
            <h1 className={styles.srOnly}>CBIcall</h1>
            <p className={styles.heroSubtitle}>
              A configuration-driven workflow dispatcher for reproducible WES,
              WGS, and mtDNA variant-calling runs.
            </p>
            <div className={styles.actions}>
              <Link className="button button--primary button--lg" to="/docs/usage/quickstart">
                Quick Start
              </Link>
              <Link className="button button--secondary button--lg" to="/docs/usage/choose-your-path">
                Choose Path
              </Link>
            </div>
          </div>
        </section>
        <section className={styles.sections}>
          <div className={styles.grid}>
            <Link className={styles.card} to="/docs/usage/choose-your-path">
              <h2>Install</h2>
              <p>Docker, Apptainer, and non-containerized setup paths.</p>
            </Link>
            <Link className={styles.card} to="/docs/usage/quickstart">
              <h2>Run</h2>
              <p>Single-sample and cohort workflows for WES, WGS, and mtDNA.</p>
            </Link>
            <Link className={styles.card} to="/docs/help/outputs">
              <h2>Inspect</h2>
              <p>Standard project layout, logs, QC files, and generated outputs.</p>
            </Link>
            <Link className={styles.card} to="/docs/technical-details/adding-a-pipeline">
              <h2>Extend</h2>
              <p>Registry-driven workflow wiring for adding new pipelines.</p>
            </Link>
          </div>
        </section>
      </main>
    </Layout>
  );
}
