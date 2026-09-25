// Optional browser verification of the script-free nested explorer.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const {pathToFileURL}=require('node:url');
const path=require('node:path'), fs=require('node:fs'), assert=require('node:assert/strict');
(async()=>{
 const html=path.resolve(process.argv[2] || 'artifacts/model_explorer/index.html');
 const browser=await chromium.launch({headless:true});
 try{
 const context=await browser.newContext({viewport:{width:1440,height:1000},offline:true,javaScriptEnabled:false});
 const page=await context.newPage(),errors=[],network=[];
 page.on('pageerror',e=>errors.push(String(e)));
 page.on('request',r=>{if(/^https?:/.test(r.url()))network.push(r.url());});
 await page.goto(pathToFileURL(html).href);
 assert.equal(await page.locator('script,input[type=search]').count(),0);
 assert.equal(await page.locator('.template-note').isVisible(),false);
 await page.screenshot({path:path.join(path.dirname(html),'overview.png'),fullPage:true});
 await page.locator('.hero a[href="#open-x"]').click();
 assert.equal(await page.locator('#node-x').getAttribute('open'),'');
 await page.locator('#node-x > summary').click();
 for(const id of ['node-x','node-x-carbon','node-x-carbon-plant','node-x-carbon-plant-tree','node-x-carbon-plant-tree-x_t_fol']){
   await page.locator('#'+id+' > summary').click();
   assert.equal(await page.locator('#'+id).getAttribute('open'),'');
 }
 const leaf=page.locator('#node-x-carbon-plant-tree-x_t_fol');
 assert.match(await leaf.locator(':scope > .inside > .equation').innerText(),/m|f_t_fol/);
 const process=leaf.locator(':scope > .inside > .nest > details').first();
 await process.locator(':scope > summary').click();
 assert.equal(await process.getAttribute('open'),'');
 assert.match(await process.locator(':scope > .inside > .equation').innerText(),/=/);
 await page.screenshot({path:path.join(path.dirname(html),'nested_leaf.png'),fullPage:true});
 await page.locator('#node-x > summary').click();
 for(const id of ['node-B','node-B-allocation','node-B-allocation-dynamic','node-B-allocation-dynamic-target','node-B-allocation-dynamic-target-target_formula']){
   await page.locator('#'+id+' > summary').click();
   assert.equal(await page.locator('#'+id).getAttribute('open'),'');
 }
 const target=page.locator('#node-B-allocation-dynamic-target-target_formula');
 assert.match(await target.locator(':scope > .inside > .equation').innerText(),/Ph−24a/);
 await target.locator(':scope > .inside > .source-detail > summary').first().click();
 assert.match(await target.locator('.source-code').first().innerText(),/f_opt_lai/);
 await page.screenshot({path:path.join(path.dirname(html),'nested_allocation.png'),fullPage:true});
 // Keyboard toggling works without JavaScript.
 await page.locator('#node-K > summary').focus();
 await page.keyboard.press('Enter');
 assert.equal(await page.locator('#node-K').getAttribute('open'),'');
 await page.setViewportSize({width:390,height:844});
 assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth),false);
 await page.screenshot({path:path.join(path.dirname(html),'mobile_nested.png'),fullPage:true});
 // Fragment navigation reveals native details ancestors in Chromium.
 await page.goto(pathToFileURL(html).href+'#eq-p_t_alloc_ass');
 assert.equal(await page.locator('#eq-p_t_alloc_ass').isVisible(),true);
 assert.deepEqual(errors,[]);assert.deepEqual(network,[]);
 const report={browser:browser.version(),javascript:false,offline:true,http_requests:network,errors,
   checks:['no search or scripts','nested x and leaf update','nested process equations','nested B and target LAI',
      'offline source','keyboard disclosure','390px overflow','definition fragment'],viewports:['1440×1000','390×844']};
 fs.writeFileSync(path.join(path.dirname(html),'browser_verification.json'),JSON.stringify(report,null,2));
 console.log(JSON.stringify(report));
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
